#!/usr/bin/env bash
# Wrapper that auto-detects a working libclang before invoking cargo.
# bindgen (used by hts-sys) reads LIBCLANG_PATH; if it points to a missing
# file the build dies. This script probes common locations and exports a
# valid path before exec'ing cargo.

set -euo pipefail

host_arch_pattern() {
    case "$(uname -m)" in
        x86_64)  echo 'x86-64' ;;
        aarch64) echo 'aarch64' ;;
        ppc64le) echo 'PowerPC' ;;
        *)       echo "$(uname -m)" ;;
    esac
}

find_libclang_in() {
    local dir="$1"
    [[ -d "$dir" ]] || return 1
    local pattern
    pattern=$(host_arch_pattern)
    local f
    for f in "$dir"/libclang.so "$dir"/libclang.so.* ; do
        [[ -e "$f" ]] || continue
        if file -L "$f" 2>/dev/null | grep -q "$pattern"; then
            echo "$dir"
            return 0
        fi
    done
    return 1
}

libclang_path_valid() {
    [[ -n "${LIBCLANG_PATH:-}" ]] && find_libclang_in "$LIBCLANG_PATH" >/dev/null
}

detect_libclang_dir() {
    local candidates=()

    if [[ -n "${CONDA_PREFIX:-}" ]]; then
        candidates+=("$CONDA_PREFIX/lib")
    fi

    if command -v llvm-config >/dev/null 2>&1; then
        local libdir
        libdir=$(llvm-config --libdir 2>/dev/null || true)
        [[ -n "$libdir" ]] && candidates+=("$libdir")
    fi

    if command -v clang >/dev/null 2>&1; then
        local clang_bin clang_prefix
        clang_bin=$(command -v clang)
        clang_prefix=$(dirname "$(dirname "$clang_bin")")
        candidates+=("$clang_prefix/lib" "$clang_prefix/lib64")
    fi

    candidates+=(
        /usr/lib/x86_64-linux-gnu
        /usr/lib64
        /usr/lib
        /usr/local/lib
    )

    local dir
    for dir in "${candidates[@]}"; do
        if find_libclang_in "$dir" >/dev/null; then
            echo "$dir"
            return 0
        fi
    done
    return 1
}

if ! libclang_path_valid; then
    if detected=$(detect_libclang_dir); then
        export LIBCLANG_PATH="$detected"
        echo "[build.sh] LIBCLANG_PATH=$LIBCLANG_PATH" >&2
    else
        echo "[build.sh] error: could not locate libclang.so in CONDA_PREFIX, llvm-config, clang prefix, or system paths" >&2
        exit 1
    fi
fi

exec cargo "$@"
