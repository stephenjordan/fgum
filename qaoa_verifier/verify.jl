#!/usr/bin/env julia
#
# verify.jl — QAOA satisfaction fraction verifier
#
# Given pre-computed QAOA angles (γ, β) and problem parameters (k, D),
# computes the exact satisfaction fraction c̃(γ, β) for QAOA on D-regular
# Max-k-XORSAT using the Basso/Villalonga branch tensor recurrence.
#
# Usage:
#   julia --project=. verify.jl                           # verify ../qaoa_data/qaoa_angles_for_verifier.csv
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

function main()
    # Default angles file: ../qaoa_data/qaoa_angles_for_verifier.csv
    default_csv = joinpath(@__DIR__, "..", "qaoa_data", "qaoa_angles_for_verifier.csv")

    if isempty(ARGS)
        if isfile(default_csv)
            return run_csv(default_csv)
        else
            print_usage()
            return 0
        end
    end

    # Check if first arg is a file
    if isfile(ARGS[1])
        return run_csv(ARGS[1])
    elseif startswith(ARGS[1], "--")
        return run_cli()
    else
        println("Error: '$(ARGS[1])' is not a file or a recognised option.")
        return 1
    end
end

function print_usage()
    println("""
QAOA Satisfaction Fraction Verifier
====================================

Usage:
  julia --project=. verify.jl                         # verify ../qaoa_data/qaoa_angles_for_verifier.csv
  julia --project=. verify.jl <angles.csv>            # verify a custom CSV
  julia --project=. verify.jl --k K --D D --gamma γ₁;γ₂;...;γₚ --beta β₁;β₂;...;βₚ [--clause-sign ±1] [--ctilde claimed]

CSV format (header line starting with # is skipped):
  k,D,p,ctilde,source,gamma,beta
  (gamma and beta fields use ; as separator for multiple values)
""")
end

function run_csv(filepath::String)
    println("╔══════════════════════════════════════════════════════════════════╗")
    println("║  QAOA Satisfaction Fraction Verifier                            ║")
    println("║  Companion to arXiv:2604.24633                                 ║")
    println("╚══════════════════════════════════════════════════════════════════╝")
    println()
    println("Reading angles from: $filepath")
    println()

    rows = parse_csv(filepath)
    if isempty(rows)
        println("No valid rows found in $filepath")
        return 1
    end

    println("Found $(length(rows)) configurations to verify.")
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
