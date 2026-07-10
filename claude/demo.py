"""Quick demo of the updated ListApparatus with ref/ambiguity columns."""
from list_apparatus import ListApparatus

app = ListApparatus.from_csv("test_collation.csv")

# ── Show the full apparatus ──
print(app)
print()

# ── Shared error groups (the clade evidence) ──
print("=" * 60)
print("SHARED ERROR GROUPS (genealogically informative)")
print("=" * 60)
for g in app.error_groups(ambiguity_levels={"certain"}):
    wits = ", ".join(sorted(g.witnesses))
    print(f"  {g.locus}: \"{g.reading}\" ≠ ref \"{g.reference}\"  →  {{{wits}}}")

# ── Singular errors (confirm independence) ──
print()
print("=" * 60)
print("SINGULAR ERRORS (individual, confirm independence)")
print("=" * 60)
for g in app.singular_errors(ambiguity_levels={"certain"}):
    wit = list(g.witnesses)[0]
    print(f"  {g.locus}: \"{g.reading}\" ≠ ref \"{g.reference}\"  →  {wit}")

# ── Frequency table for Froger ──
print()
print("=" * 60)
print("GROUP FREQUENCIES (input for Froger)")
print("=" * 60)
freq = app.group_by_frequency(
    only_errors=True,
    min_group_size=2,
    ambiguity_levels={"certain"},
)
for witnesses, data in sorted(freq.items(), key=lambda x: len(x[0])):
    wits = ", ".join(sorted(witnesses))
    loci = ", ".join(data["loci"])
    print(f"  {{{wits}}}  ×{data['frequency']}  (loci: {loci})")

# ── Laminar check ──
print()
print("=" * 60)
print("LAMINAR FAMILY CHECK")
print("=" * 60)
is_laminar, conflicts = app.check_laminar(ambiguity_levels={"certain"})
if is_laminar:
    print("  ✓ All shared error groups are pairwise compatible.")
    print("  → No contamination detected. Clean tree can be built.")
else:
    print(f"  ✗ {len(conflicts)} incompatible pair(s) found:")
    for a, b in conflicts:
        print(f"    {set(a)}  ⟷  {set(b)}")