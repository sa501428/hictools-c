# Future TODO

These are deferred edge cases and interoperability hardening items from the
final V10 writer and V9-to-V10 converter review. They are not blockers for the
ordinary count-based conversion path.

## Score-matrix type inference during conversion

Infer value types for an entire derived source/target family before writing it.
A fractional score source can aggregate exactly to an integral-looking target;
independent per-resolution inference currently rejects that valid conversion
unless the user supplies `--scores`.

Add coverage for both 2x and 5x derived score resolutions, including exact
float-bit comparison of the preserved target.

## Empty matrices with forced score semantics

Clarify whether `--scores` declares score semantics globally or only for
chromosome pairs present in the input. If it is global, emit explicit empty
`SCORE_FLOAT32` matrix records for otherwise missing canonical chromosome pairs,
because an omitted V10 pair is defined as an empty `COUNT_UINT` matrix.

## Aggregation.NONE interoperability

Allow the internal V10 reader used by `addnorm` to accept `Aggregation.NONE` on
materialized resolutions. Continue to require `Aggregation.SUM` for every
derived resolution. Files written by hictools-c currently use `SUM`, so this is
primarily compatibility with other conforming V10 writers.

## Defensive conformance validation

- Validate UTF-8 for every V10 `cstr` written or read.
- Validate strict resolution ordering in the internal reader.
- Validate that indexed matrix-block intervals do not overlap and require a
  stored block length greater than the 16-byte `H10B` header.
- Fully validate caller-supplied derived FRAG records in the generic writer.
