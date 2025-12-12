#!/usr/bin/env julia

using Base.Threads
using Printf
using Statistics
using ArgParse
using MurmurHash3

const KMER_HASH_SEED = UInt32(460)


# ============ argument parsing ============ #

mutable struct Options
    indir::String
    outdir::String
    kmer::String
    taxa::String
    number::Int
    procs::Int
    canonical::Int
    kmerlim::Float64
    presabs::Int
    hashing::Int
end


function get_options()::Options
    s = ArgParseSettings()

    @add_arg_table s begin
        "-i"
            help = "Input directory"
            arg_type = String
            required = true
        "-o"
            help = "Output directory"
            arg_type = String
            required = true
        "-k"
            help = "K-mer size or range (e.g., 5 or 1-50)"
            arg_type = String
            required = true
        "-t"
            help = "Taxon identifier"
            arg_type = String
            required = true
        "-n"
            help = "Sampling rate"
            arg_type = Int
            required = true
        "-p"
            help = "Number of processes/threads"
            arg_type = Int
            default = 1
        "-c"
            help = "Canonical kmers (0/1)"
            arg_type = Int
            default = 0
        "-s"
            help = "K-mer limit"
            arg_type = Float64
            default = 1e99
        "-a"
            help = "Presence/absence mode (0/1)"
            arg_type = Int
            default = 0
        "-x"
            help = "Hashing mode (1=hash, 2=raw)"
            arg_type = Int
            default = 1
    end

    parsed = ArgParse.parse_args(ARGS, s)

    return Options(
        parsed["i"],
        parsed["o"],
        parsed["k"],
        parsed["t"],
        parsed["n"],
        parsed["p"],
        parsed["c"],
        parsed["s"],
        parsed["a"],
        parsed["x"],
    )
end

# ============ FASTA reading (handles gz via external gzip -dc) ============ #

function open_text_maybe_gzip(path::String)
    if endswith(path, ".gz")
        # pipe from gzip -dc
        p = open(`gzip -dc $path`)
        return p
    else
        return open(path, "r")
    end
end

#Returns all concatenated sequences (uppercase) from a FASTA file.
function parse_fasta_sequences(path::String)::Vector{String}
    io = open_text_maybe_gzip(path)
    sequences = String[]
    current = IOBuffer()

    firstseq = true
    for line in eachline(io)
        line = strip(line)
        isempty(line) && continue
        if startswith(line, '>')
            # header line
            if position(current) > 0
                push!(sequences, uppercase(String(take!(current))))
            end
        else
            print(current, line)
        end
    end
    if position(current) > 0
        push!(sequences, uppercase(String(take!(current))))
    end
    close(io)
    return sequences
end

# ============ DNA helper functions ============ #

const AMBIG_RE = r"[BDEFHIJKLMNOPQRSUVWXYZ]"

function revdnacomp(dna::AbstractString)::String
    trans = Dict(
        'A'=>'T','C'=>'G','G'=>'C','T'=>'A',
        'a'=>'t','c'=>'g','g'=>'c','t'=>'a'
    )
    buf = IOBuffer()
    for c in reverse(dna)
        print(buf, get(trans, c, c))
    end
    return String(take!(buf))
end

# ============ Entropy ============ #

"""
    ent(dat::Vector{Float64}) -> (entropy, normalized_entropy)

Uses natural log, matching the Perl ent() sub.
"""
function ent(dat::Vector{Float64})
    entropy = 0.0
    c = 0
    for d in dat
        d == 0.0 && continue
        c += 1
        entropy += d * log(d)
    end
    entval = -entropy
    nent = c > 0 ? -entropy / log(c) : 0.0
    return entval, nent
end

# ============ KL divergence ============ #

"""
    kl_divergence(g, kmers, kpg, genomes, allk)

kmers :: Dict{String, Dict{Int, Int}}
kpg   :: Dict{Int, Int}
genomes :: Vector{Int} (indices)
"""
function kl_divergence(
    g::Int,
    kmers::Dict{String, Dict{Int, Int}},
    kpg::Dict{Int, Int},
    genomes::Vector{Int},
    allk::Int
)::Float64
    klval = 0.0
    kpg_g = kpg[g]
    for (_, counts_by_gen) in kmers
        psi = get(counts_by_gen, g, 0) / kpg_g
        psi == 0.0 && continue
        pts = 0
        for subg in genomes
            pts += get(counts_by_gen, subg, 0)
        end
        pss = pts / allk
        kli = psi * log(psi / pss)
        klval += kli
    end
    return klval
end

# ============ Hashing ============ #

# MurmurHash3 x64_128 seeded hash, returns low 64 bits.
function hash_kmer(cfrag::AbstractString)::UInt64
    h1, h2 = mmhash128_c(cfrag, KMER_HASH_SEED)
    return UInt64(h1)  # use low 64 bits for sampling and output
end


# ============ Genome k-mer hashing ============ #

struct KmerJob
    indir::String
    outdir_k::String
    genome::String
    kmersize::Int
    canonical::Int
    hashing::Int
    number::Int
    kmer_limit::Float64
    thresh::Float64
    motion::Int
end

function process_genome_to_hashes(job::KmerJob)
    indir, outdir_k, genome, kmersize, canonical, hashing, number, kmer_limit, thresh, motion =
        job.indir, job.outdir_k, job.genome, job.kmersize, job.canonical,
        job.hashing, job.number, job.kmer_limit, job.thresh, job.motion

    gzpath = joinpath(indir, genome * ".gz")
    fasta_path = isfile(gzpath) ? gzpath : joinpath(indir, genome)

    out_path = joinpath(outdir_k, genome * ".hashes")
    kc = 0
    signal = false

    seqs = parse_fasta_sequences(fasta_path)
    open(out_path, "w") do out_io
        for seq in seqs
            signal && break
            seqlen = lastindex(seq)
            idx = firstindex(seq)
            while idx + kmersize - 1 <= seqlen
                ffrag = @views seq[idx : idx + kmersize - 1]

                if occursin(AMBIG_RE, ffrag)
                    idx += motion
                    continue
                end

                # canonicalization
                cfrag = if canonical == 1
                    rc = revdnacomp(ffrag)
                    ffrag < rc ? String(ffrag) : rc
                else
                    String(ffrag)
                end

                if hashing == 1 && kmersize > 3
                    h = hash_kmer(cfrag)
                    if h < thresh
                        println(out_io, h)
                        kc += 1
                    end
                elseif hashing == 1 && kmersize <= 3
                    # subsample small k-mers
                    if rand() <= 1.0 / number
                        println(out_io, cfrag)
                        kc += 1
                    end
                else
                    # no hashing: keep all kmers
                    println(out_io, cfrag)
                    kc += 1
                end

                idx += motion
                if kc >= kmer_limit
                    signal = true
                    break
                end
            end
        end
    end
    return nothing
end

# ============ Genome length computation ============ #

function compute_genome_lengths(indir::String, genomes::Vector{String})
    lengths = Dict{String, Int}()
    for g in genomes
        gzpath = joinpath(indir, g * ".gz")
        fasta_path = isfile(gzpath) ? gzpath : joinpath(indir, g)
        total_len = 0
        for seq in parse_fasta_sequences(fasta_path)
            total_len += ncodeunits(seq)
        end
        lengths[g] = total_len
    end
    return lengths
end

# ============ Main ============ #

function main()
    opt = get_options()
    println("Using $(nthreads()) Julia threads.")

    indir     = opt.indir
    outdir    = opt.outdir
    kmer_arg  = opt.kmer
    taxa      = opt.taxa
    number    = opt.number
    canonical = opt.canonical
    kmerlim   = opt.kmerlim
    presabs   = opt.presabs
    hashing   = opt.hashing

    isdir(outdir) || mkpath(outdir)

    # Discover genomes
    genomes = String[]
    for entry in readdir(indir)
        startswith(entry, ".") && continue
        base = endswith(entry, ".gz") ? entry[1:end-3] : entry
        if !(base in genomes)
            push!(genomes, base)
        end
    end
    sort!(genomes)

    # Index genomes as Int for faster dicts
    genome_index = Dict{String, Int}()
    for (i, g) in enumerate(genomes)
        genome_index[g] = i
    end

    lengths = compute_genome_lengths(indir, genomes)

    maxhash = 2.0^64
    thresh = maxhash / number
    motion = 1

    ents  = Dict{Int, Float64}()
    nents = Dict{Int, Float64}()
    effg  = Dict{Int, Float64}()
    effk  = Dict{Int, Float64}()

    # Parse kmer argument
    ksizes = Int[]
    if occursin("-", kmer_arg)
        parts = split(kmer_arg, "-")
        s = parse(Int, parts[1])
        e = parse(Int, parts[2])
        append!(ksizes, s:e)
    else
        push!(ksizes, parse(Int, kmer_arg))
    end
    sort!(ksizes)

    # === loop over k sizes === #
    for kmersize in ksizes
        outdir_k = joinpath(outdir, string(kmersize))
        isdir(outdir_k) || mkpath(outdir_k)

        # k-mer counting in parallel
        jobs = KmerJob[]
        for g in genomes
            push!(jobs, KmerJob(indir, outdir_k, g, kmersize,
                                 canonical, hashing, number, kmerlim, thresh, motion))
        end

        @info "Counting kmers at k=$kmersize for $(length(genomes)) genomes"
        # Threaded map
        @threads for j in 1:length(jobs)
            job = jobs[j]
            @info "  Counting $(job.genome)"
            process_genome_to_hashes(job)
        end

        # Collate hashes
        kmers = Dict{String, Dict{Int, Int}}()
        kmersperg = Dict{Int, Int}()
        allkmers = 0

        mast_path = joinpath(outdir_k, "stats.out")
        open(mast_path, "w") do mast
            for g in genomes
                gi = genome_index[g]
                @info "Collating hashes for $g (k=$kmersize)"
                gkmerscount = 0
                gkmers = Dict{String, Int}()
                cache = Set{String}()

                hashes_path = joinpath(outdir_k, g * ".hashes")
                isfile(hashes_path) || continue

                open(hashes_path, "r") do hio
                    if presabs == 1
                        for line in eachline(hio)
                            hv = strip(line)
                            isempty(hv) && continue
                            if !(hv in cache)
                                push!(cache, hv)
                                gkmerscount += 1
                                gkmers[hv] = get(gkmers, hv, 0) + 1
                                inner = get!(kmers, hv, Dict{Int, Int}())
                                inner[gi] = get(inner, gi, 0) + 1
                                kmersperg[gi] = get(kmersperg, gi, 0) + 1
                                allkmers += 1
                            end
                        end
                    else
                        for line in eachline(hio)
                            hv = strip(line)
                            isempty(hv) && continue
                            gkmerscount += 1
                            gkmers[hv] = get(gkmers, hv, 0) + 1
                            inner = get!(kmers, hv, Dict{Int, Int}())
                            inner[gi] = get(inner, gi, 0) + 1
                            kmersperg[gi] = get(kmersperg, gi, 0) + 1
                            allkmers += 1
                        end
                    end
                end

                if gkmerscount > 0
                    gkfreqs = Float64[]
                    for cnt in values(gkmers)
                        push!(gkfreqs, cnt / gkmerscount)
                    end
                    gkent, gknent = ent(gkfreqs)
                    egkent = exp(gkent)
                    uniqgkmers = length(gkmers)
                    @printf(mast, "%s\t%d\t%d\t%.6f\t%.6f\t%.6f\n",
                            g, gkmerscount, uniqgkmers, gknent, gkent, egkent)
                end
            end

            # Ensemble entropies
            allkfreqs = Float64[]
            for inner in values(kmers)
                kcount = 0
                for v in values(inner)
                    kcount += v
                end
                push!(allkfreqs, kcount / allkmers)
            end
            kent, knent = ent(allkfreqs)
            ekent = exp(kent)
            uniqkmers = length(kmers)

            ents[kmersize]  = kent
            nents[kmersize] = knent
            effk[kmersize]  = ekent

            @printf(mast, "ENSEMBLE\t%d\t%d\t%.6f\t%.6f\t%.6f\n",
                    allkmers, uniqkmers, knent, kent, ekent)

            # KL divergence + beta entropy
            genomes_idx = collect(values(genome_index))  # 1..N
            betaent = 0.0

            # Parallel KL in threads (use a Vector for thread-safe writes)
            max_idx = maximum(genomes_idx)
            bevals = zeros(Float64, max_idx)

            @threads for gidx in genomes_idx
                if !haskey(kmersperg, gidx)
                    continue
                end
                gname = genomes[gidx]
                @info "Calculating KL for $gname (k=$kmersize)"
                klv = kl_divergence(gidx, kmers, kmersperg, genomes_idx, allkmers)
                weight = kmersperg[gidx] / allkmers
                be = klv * weight
                bevals[gidx] = be  # each index is written by only one thread

                be_path = joinpath(outdir_k, gname * ".betaent")
                open(be_path, "w") do beio
                    @printf(beio, "%s\t%.6f\t%.6f\t%.6f\n", gname, be, klv, weight)
                end
            end

            for v in bevals
                betaent += v
            end
            betadiv = exp(betaent)
            effg[kmersize] = betadiv

            @printf(mast, "betaEnt=%.6f\n", betaent)
            @printf(mast, "betaDiv=%.6f\n", betadiv)
        end

        # remove hash files
        for g in genomes
            hashes_path = joinpath(outdir_k, g * ".hashes")
            isfile(hashes_path) && rm(hashes_path)
        end
    end

    # === Post-processing: curves & structure === #
    maxksize = maximum(ksizes)

    # delta H
    deltaents = Dict{Int, Float64}()
    entcache = 0.0
    for k in ksizes
        d = ents[k] - entcache
        deltaents[k] = d
        entcache = ents[k]
    end

    exent = 0.0
    pred = 0.0
    exentstep = Dict{Int, Float64}()
    deltadeltaents = Dict{Int, Float64}()
    deltaentcache = 0.0
    for k in ksizes
        d2 = deltaents[k] - deltaentcache
        deltadeltaents[k] = d2
        exent += (deltaents[k] - deltaents[maxksize])
        exentstep[k] = exent
        pred += d2
        deltaentcache = deltaents[k]
    end

    transinfo = 0.0
    transinfostep = Dict{Int, Float64}()
    for k in ksizes
        transinfo += (exent + deltaents[maxksize] * k) - ents[k]
        transinfostep[k] = transinfo
    end

    ent_hills_path = joinpath(outdir, "ent_hills.out")
    open(ent_hills_path, "w") do out
        header = join([
            "kmersize",
            outdir * "_H",
            outdir * "_nH",
            outdir * "_dH",
            outdir * "d2H",
            outdir * "exentStep",
            outdir * "transinfoStep",
            outdir * "_effKmers",
            outdir * "_khill",
        ], '\t')
        println(out, header)

        for k in ksizes
            @printf(out, "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    k,
                    ents[k],
                    nents[k],
                    deltaents[k],
                    deltadeltaents[k],
                    exentstep[k],
                    transinfostep[k],
                    effk[k],
                    effg[k])
        end
    end

    # structure.out
    lenvals = collect(values(lengths))
    totlen = sum(lenvals)
    meanlen = isempty(lenvals) ? 0.0 : mean(lenvals)
    stdv = length(lenvals) > 1 ? std(lenvals) : 0.0

    maxh = isempty(ents) ? 0.0 : maximum(values(ents))
    maxeffk = isempty(effk) ? 0.0 : maximum(values(effk))
    maxeffg = isempty(effg) ? 0.0 : maximum(values(effg))

    minnh = Inf
    minnhk = 0
    mind2h = Inf
    mind2hk = 0

    for k in ksizes
        nh = get(nents, k, Inf)
        d2 = get(deltadeltaents, k, Inf)
        if nh < minnh
            minnh = nh
            minnhk = k
        end
        if d2 < mind2h
            mind2h = d2
            mind2hk = k
        end
    end
    minnh = isfinite(minnh) ? minnh : 0.0
    mind2h = isfinite(mind2h) ? mind2h : 0.0

    structure_path = joinpath(outdir, "structure.out")
    open(structure_path, "w") do st
        @printf(st,
                "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%d\n",
                taxa,
                totlen,
                meanlen,
                stdv,
                deltaents[maxksize],
                exent,
                transinfo,
                pred,
                maxh,
                maxeffk,
                maxeffg,
                minnh,
                minnhk,
                mind2h,
                mind2hk)
    end
end

main()
