#!/bin/bash

# Configuration: Use LMP_EXE from environment or default to local build path
LMP_EXE=${LMP_EXE:-../../../build-grace/lmp}

echo "Using LAMMPS executable: $LMP_EXE"

# Function to run a command and check for failure
check_run() {
    local cmd=$1
    local screen_log=$2
    eval "$cmd" > "$screen_log" 2>&1
    if [ $? -ne 0 ]; then
        cat "$screen_log"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "           RUN ERROR"
        echo "   COMMAND FAILED: $cmd"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        exit 1
    fi
}

# Function to extract and normalize thermo output for comparison
extract_thermo() {
    local logfile=$1
    # Extract Step header and the following data line (ignoring trailing whitespace)
    grep -A 1 "Step[[:space:]]\+TotEng" "$logfile" | head -n 2 | sed 's/[[:space:]]\+$//'
}

# Default to fast mode (chunked only)
FAST=1
if [[ "$*" == *"--full"* ]]; then
    FAST=0
fi

SERIAL_LOG="log-serial.lammps"
SERIAL_SCREEN="screen-serial.log"
CHUNKED_LOG="log-chunked.lammps"
CHUNKED_SCREEN="screen-chunked.log"
MPI_LOG="log-mpi.lammps"
MPI_SCREEN="screen-mpi.log"

if [ $FAST -eq 0 ]; then
    echo "Running Serial simulation (grace)..."
    check_run "$LMP_EXE -in in-bcc-serial.lammps -log $SERIAL_LOG" "$SERIAL_SCREEN"
    echo "Done."
else
    echo "Skipping Serial simulation (--fast is default, use --full to include it)."
fi

echo "Running Chunked simulation (grace/1layer/chunk)..."
check_run "$LMP_EXE -in in-bcc-chunked.lammps -log $CHUNKED_LOG" "$CHUNKED_SCREEN"
mv structure_chunked.dump structure_chunked_serial.dump
echo "Done."

echo "Running MPI Chunked simulation (grace/1layer/chunk -np 2)..."
check_run "mpirun -np 2 $LMP_EXE -in in-bcc-chunked.lammps -log $MPI_LOG" "$MPI_SCREEN"
mv structure_chunked.dump structure_chunked_mpi.dump
echo "Done."

extract_thermo "$SERIAL_LOG" > thermo-serial.tmp
echo "--- Thermo Output: Serial ---"
cat thermo-serial.tmp
grep "Loop time of" "$SERIAL_LOG" | sed 's/[[:space:]]\+$//'
echo "------------------------------"

extract_thermo "$CHUNKED_LOG" > thermo-chunked.tmp
echo "--- Thermo Output: Chunked ---"
cat thermo-chunked.tmp
grep "Loop time of" "$CHUNKED_LOG" | sed 's/[[:space:]]\+$//'
echo "-------------------------------"

extract_thermo "$MPI_LOG" > thermo-mpi.tmp
echo "--- Thermo Output: MPI ---"
cat thermo-mpi.tmp
grep "Loop time of" "$MPI_LOG" | sed 's/[[:space:]]\+$//'
echo "-------------------------------"

# Compare thermo
THERMO_MATCH=true
if ! diff thermo-serial.tmp thermo-chunked.tmp > /dev/null; then
    THERMO_MATCH=false
    echo "Thermodynamic data mismatch between Serial and Chunked!"
    diff -u thermo-serial.tmp thermo-chunked.tmp
fi

if ! diff thermo-serial.tmp thermo-mpi.tmp > /dev/null; then
    THERMO_MATCH=false
    echo "Thermodynamic data mismatch between Serial and MPI!"
    diff -u thermo-serial.tmp thermo-mpi.tmp
fi

[ "$THERMO_MATCH" = true ] && echo "Thermodynamic data matches across all runs."

# Compare dumps
DUMP_MATCH=true
if ! diff structure_serial.dump structure_chunked_serial.dump > /dev/null; then
    DUMP_MATCH=false
    echo "Dump files mismatch between Serial and Chunked Serial!"
    diff -u structure_serial.dump structure_chunked_serial.dump | head -n 20
fi

if ! diff structure_serial.dump structure_chunked_mpi.dump > /dev/null; then
    DUMP_MATCH=false
    echo "Dump files mismatch between Serial and MPI!"
    diff -u structure_serial.dump structure_chunked_mpi.dump | head -n 20
fi

[ "$DUMP_MATCH" = true ] && echo "Dump files match across all runs."

# Cleanup
rm thermo-serial.tmp thermo-chunked.tmp thermo-mpi.tmp

# Final result
echo ""
if [ "$THERMO_MATCH" = true ] && [ "$DUMP_MATCH" = true ]; then
    echo "================================"
    echo "   ALL TESTS SUCCESSFUL"
    echo "================================"
    exit 0
else
    echo "================================"
    echo "   TEST FAILED: MISMATCH"
    echo "================================"
    exit 1
fi
