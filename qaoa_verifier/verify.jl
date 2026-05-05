#!/usr/bin/env julia
#
# verify.jl — QAOA satisfaction fraction verifier
#
# Given pre-computed QAOA angles (γ, β) and problem parameters (k, D),
# computes the exact satisfaction fraction c̃(γ, β) for QAOA on D-regular
# Max-k-XORSAT using the Basso/Villalonga branch tensor recurrence.
#
# Usage:
#   julia --project=. verify.jl                           # verify all 254 configurations
#   julia --project=. verify.jl --p 10                    # verify only depth p=10 (16 configs)
#   julia --project=. verify.jl --max-p 12                # verify depths p=1..12
#   julia --project=. verify.jl --filter-k 3 --filter-D 4 # verify only (k=3, D=4) across all depths
#   julia --project=. verify.jl path/to/angles.csv        # verify a custom CSV
#   julia --project=. verify.jl --k 3 --D 4 --gamma "0.55;0.72" --beta "0.36;0.21"
#
# CSV format (header/comment lines starting with # are skipped):
#   k,D,p,ctilde,source,gamma,beta
#   (gamma and beta fields use ; as separator for multiple depth values)
#
# Clause sign is inferred from k: k=2 → clause_sign=-1 (MaxCut), k≥3 → clause_sign=+1 (XORSAT).
#
# Author: John S Azariah — Centre for Quantum Software and Information, UTS

using QaoaVerifier
import Printf

# ──────────────────────────────────────────────────────────────────────────────
# CSV parsing
# ──────────────────────────────────────────────────────────────────────────────

function parse_angles(field::AbstractString)::Vector{Float64}
    parse.(Float64, split(strip(field), ';'))
end

function infer_clause_sign(k::Int)
    k == 2 ? -1 : 1
end

struct VerificationRow
    k::Int
    D::Int
    p::Int
    ctilde_claimed::Float64
    gamma::Vector{Float64}
    beta::Vector{Float64}
    clause_sign::Int
    source::String
end

function parse_csv(filepath::String)::Vector{VerificationRow}
    rows = VerificationRow[]
    for line in eachline(filepath)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, '#') && continue

        fields = split(stripped, ',')
        length(fields) ≥ 7 || continue

        k = parse(Int, fields[1])
        D = parse(Int, fields[2])
        p = parse(Int, fields[3])
        ctilde = parse(Float64, fields[4])
        source = String(fields[5])
        gamma = parse_angles(fields[6])
        beta = parse_angles(fields[7])

        length(gamma) == p || error("Row k=$k D=$D p=$p: gamma has $(length(gamma)) entries, expected $p")
        length(beta) == p || error("Row k=$k D=$D p=$p: beta has $(length(beta)) entries, expected $p")

        cs = infer_clause_sign(k)
        push!(rows, VerificationRow(k, D, p, ctilde, gamma, beta, cs, source))
    end
    rows
end

# ──────────────────────────────────────────────────────────────────────────────
# CLI argument parsing
# ──────────────────────────────────────────────────────────────────────────────

function parse_cli_args(args)
    k = nothing; D = nothing
    gamma = nothing; beta = nothing
    clause_sign = nothing
    ctilde_claimed = nothing

    i = 1
    while i ≤ length(args)
        arg = args[i]
        if arg == "--k"
            k = parse(Int, args[i+1]); i += 2
        elseif arg == "--D"
            D = parse(Int, args[i+1]); i += 2
        elseif arg == "--gamma"
            gamma = parse_angles(args[i+1]); i += 2
        elseif arg == "--beta"
            beta = parse_angles(args[i+1]); i += 2
        elseif arg == "--clause-sign"
            clause_sign = parse(Int, args[i+1]); i += 2
        elseif arg == "--ctilde"
            ctilde_claimed = parse(Float64, args[i+1]); i += 2
        else
            i += 1
        end
    end

    (k=k, D=D, gamma=gamma, beta=beta, clause_sign=clause_sign, ctilde_claimed=ctilde_claimed)
end

# ──────────────────────────────────────────────────────────────────────────────
# Formatting
# ──────────────────────────────────────────────────────────────────────────────

function problem_label(k::Int, clause_sign::Int)
    if k == 2 && clause_sign == -1
        return "MaxCut"
    elseif clause_sign == 1
        return "Max-$(k)-XORSAT"
    else
        return "k=$(k) cs=$(clause_sign)"
    end
end

function format_result(k, D, p, ctilde, ctilde_claimed, elapsed, clause_sign)
    label = problem_label(k, clause_sign)
    status = if isnan(ctilde)
        "FAIL (NaN)"
    elseif ctilde_claimed !== nothing
        diff = abs(ctilde - ctilde_claimed)
        if diff < 1e-6
            "OK (Δ=$(Printf.@sprintf("%.2e", diff)))"
        else
            "MISMATCH (Δ=$(Printf.@sprintf("%.2e", diff)))"
        end
    else
        ""
    end

    claimed_str = ctilde_claimed !== nothing ? Printf.@sprintf("  claimed=%.12f", ctilde_claimed) : ""
    Printf.@sprintf("  (k=%d, D=%d, p=%2d) [%-16s]  c̃=%.12f%s  %s  [%.1fs]",
        k, D, p, label, ctilde, claimed_str, status, elapsed)
end

# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

function parse_filter_args(args)
    filter_p = nothing
    max_p = nothing
    filter_k = nothing
    filter_D = nothing
    remaining = String[]

    i = 1
    while i ≤ length(args)
        arg = args[i]
        if arg == "--p"
            filter_p = parse(Int, args[i+1]); i += 2
        elseif arg == "--max-p"
            max_p = parse(Int, args[i+1]); i += 2
        elseif arg == "--filter-k"
            filter_k = parse(Int, args[i+1]); i += 2
        elseif arg == "--filter-D"
            filter_D = parse(Int, args[i+1]); i += 2
        else
            push!(remaining, arg); i += 1
        end
    end

    (filter_p=filter_p, max_p=max_p, filter_k=filter_k, filter_D=filter_D, remaining=remaining)
end

function should_include(row, filters)
    filters.filter_p !== nothing && row.p != filters.filter_p && return false
    filters.max_p !== nothing && row.p > filters.max_p && return false
    filters.filter_k !== nothing && row.k != filters.filter_k && return false
    filters.filter_D !== nothing && row.D != filters.filter_D && return false
    true
end

function main()
    # Default angles file: ../qaoa_data/qaoa_angles_for_verifier.csv
    default_csv = joinpath(@__DIR__, "..", "qaoa_data", "qaoa_angles_for_verifier.csv")

    filters = parse_filter_args(ARGS)
    has_filters = filters.filter_p !== nothing || filters.max_p !== nothing ||
                  filters.filter_k !== nothing || filters.filter_D !== nothing
    remaining = filters.remaining

    if isempty(remaining)
        if isfile(default_csv)
            return run_csv(default_csv; filters)
        elseif !has_filters
            print_usage()
            return 0
        else
            println("Error: no angles CSV found at $default_csv")
            return 1
        end
    end

    # Check if first arg is a file
    if isfile(remaining[1])
        return run_csv(remaining[1]; filters)
    elseif startswith(remaining[1], "--")
        return run_cli()
    else
        println("Error: '$(remaining[1])' is not a file or a recognised option.")
        return 1
    end
end

function print_usage()
    println("""
QAOA Satisfaction Fraction Verifier
====================================

Usage:
  julia --project=. verify.jl                         # verify all 254 configurations
  julia --project=. verify.jl --p 10                   # verify only depth p=10
  julia --project=. verify.jl --max-p 12               # verify depths p=1..12
  julia --project=. verify.jl --filter-k 3 --filter-D 4  # verify only (k=3, D=4)
  julia --project=. verify.jl <angles.csv>             # verify a custom CSV
  julia --project=. verify.jl --k K --D D --gamma γ₁;γ₂;...;γₚ --beta β₁;β₂;...;βₚ [--clause-sign ±1] [--ctilde claimed]

CSV format (header line starting with # is skipped):
  k,D,p,ctilde,source,gamma,beta
  (gamma and beta fields use ; as separator for multiple values)
""")
end

function run_csv(filepath::String; filters=(
        filter_p=nothing, max_p=nothing, filter_k=nothing, filter_D=nothing,
        remaining=String[]))
    println("╔══════════════════════════════════════════════════════════════════╗")
    println("║  QAOA Satisfaction Fraction Verifier                            ║")
    println("║  Companion to arXiv:2604.24633                                 ║")
    println("╚══════════════════════════════════════════════════════════════════╝")
    println()
    println("Reading angles from: $filepath")

    all_rows = parse_csv(filepath)
    if isempty(all_rows)
        println("No valid rows found in $filepath")
        return 1
    end

    rows = filter(r -> should_include(r, filters), all_rows)
    if isempty(rows)
        println("No configurations match the given filters.")
        return 1
    end

    filter_desc = String[]
    filters.filter_p !== nothing && push!(filter_desc, "p=$(filters.filter_p)")
    filters.max_p !== nothing && push!(filter_desc, "p≤$(filters.max_p)")
    filters.filter_k !== nothing && push!(filter_desc, "k=$(filters.filter_k)")
    filters.filter_D !== nothing && push!(filter_desc, "D=$(filters.filter_D)")
    if !isempty(filter_desc)
        println("Filter: $(join(filter_desc, ", "))")
    end
    println()
    println("Verifying $(length(rows)) of $(length(all_rows)) configurations.")
    println("─" ^ 100)

    any_failed = false
    any_mismatch = false
    for row in rows
        elapsed = @elapsed ctilde = verify_satisfaction_fraction(
            row.k, row.D, row.gamma, row.beta;
            clause_sign = row.clause_sign,
        )

        result_str = format_result(row.k, row.D, row.p, ctilde, row.ctilde_claimed, elapsed, row.clause_sign)
        println(result_str)

        if isnan(ctilde) || !isfinite(ctilde)
            any_failed = true
        end
        if abs(ctilde - row.ctilde_claimed) ≥ 1e-6
            any_mismatch = true
        end
    end

    println("─" ^ 100)
    if any_failed
        println("⚠  Some evaluations produced NaN/Inf — check parameters.")
        return 1
    elseif any_mismatch
        println("⚠  Some evaluations did not match claimed values (Δ ≥ 1e-6).")
        return 1
    else
        println("✓  All $(length(rows)) evaluations verified successfully.")
        return 0
    end
end

function run_cli()
    cli = parse_cli_args(ARGS)
    if cli.k === nothing || cli.D === nothing || cli.gamma === nothing || cli.beta === nothing
        println("Error: --k, --D, --gamma, and --beta are all required in CLI mode.")
        return 1
    end

    cs = cli.clause_sign !== nothing ? cli.clause_sign : infer_clause_sign(cli.k)
    p = length(cli.gamma)

    println("╔══════════════════════════════════════════════════════════════════╗")
    println("║  QAOA Satisfaction Fraction Verifier                            ║")
    println("╚══════════════════════════════════════════════════════════════════╝")
    println()

    elapsed = @elapsed ctilde = verify_satisfaction_fraction(
        cli.k, cli.D, cli.gamma, cli.beta;
        clause_sign = cs,
    )

    println(format_result(cli.k, cli.D, p, ctilde, cli.ctilde_claimed, elapsed, cs))

    if isnan(ctilde) || !isfinite(ctilde)
        return 1
    end
    return 0
end

exit(main())
