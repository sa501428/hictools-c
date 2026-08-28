# Hi-C binary short (HBS), version 1

HBS is a streaming, gzip-compressed interchange format for raw integer contacts
at one BP resolution. Its extension is `.hbs.gz`. It is independent of the
V9/V10 hic wire formats and of straw's older HICSLICE format.

All multi-byte integers are **unsigned little-endian**. Fields are serialized
individually; there is no alignment padding. The layout below describes the
**decompressed** byte stream. Gzip is required, including its integrity trailer.

## Header

| Field | Bytes | Value |
|---|---:|---|
| Magic | 8 | Hex `48 49 43 42 53 00 0d 0a` (`HICBS\0\r\n`) |
| Version | 2 | `1` |
| Reserved flags | 2 | `0`; other values are rejected |
| BP resolution | 4 | Positive, at most `INT32_MAX` |
| Number of chromosomes | 4 | 0–65,536 |

Each chromosome table entry contains:

| Field | Bytes | Meaning |
|---|---:|---|
| Name length | 2 | 1–4096 bytes |
| Name | name length | Name bytes, without a NUL terminator or embedded NUL |
| Chromosome length | 8 | Positive BP length |

IDs are zero-based table indices. Names must be unique. The table excludes the
synthetic ALL overview; no particular real-chromosome ordering is required.
The complete header is limited to 16 MiB. Zero chromosomes are allowed for an
empty export.

## Contact records

| Field | Bytes |
|---|---:|
| Chromosome 1 ID | 2 |
| Bin 1 index | 4 |
| Chromosome 2 ID | 2 |
| Bin 2 index | 4 |
| Count or escape | 2 |

A count from 0 through 65,534 is stored directly. A count field of **65,535**
is followed by an **8-byte uint64 count**, which must be at least 65,535.
Thus ordinary records use **14 bytes**, and escaped records use **22 bytes**.
Exporters omit zero counts; readers accept them.

Positions are zero-based bin starts: `position = bin_index * resolution`.
Each start must be strictly less than its chromosome length. The final bin may
extend past that length. Counts are raw observed NONE integer contacts, not
normalization weights, O/E values, fragment counts, or floating-point scores.

There is no record count or terminator: records continue to the end of the
decompressed stream. A partial record, incomplete escaped count, truncated gzip
stream, or checksum failure is an error.

For preprocessing, each unordered chromosome pair must occupy one contiguous
block of records. Records within a pair need not be sorted. Straw emits pairs
in source chromosome order and canonicalizes cis coordinates, without expanding
the reflected matrix triangle. Repeated cells can be aggregated by the builder.

## Producer and consumer behavior

Straw writes through a temporary file beside the destination, closes gzip, and
renames only on success. It rejects input/output aliases. Unsupported raw
fractional, negative, or nonfinite values fail instead of being rounded.

hictools-c detects `.hbs.gz` automatically, or accepts `-f hbs`. It maps names
through the supplied genome (including its existing chromosome aliases), checks
lengths, and rejects unknown chromosomes or multiple table entries mapping to
the same genome chromosome. The genome may have extra chromosomes.

Rebuild resolutions must be positive multiples of the source resolution.
The current builders require reconstructed BP positions to fit signed int32;
the wire format itself can represent larger positions. V9 uses its existing
float32 weight/storage path. V10 preserves HBS uint64 counts through parsing,
spooling, and checked aggregation; explicit `--scores` requests float32 storage.
Expected/normalization vectors remain floating point.

