# Feature coverage matrix

Single source of truth for **what works** in survatr. Mirrors causatr's
`FEATURE_COVERAGE_MATRIX.md` convention: every PR that adds, removes, or
changes a feature MUST update this file and the corresponding tests.

## Legend

- 🟢 — supported, truth-based test pinned against an analytical or external
  reference (`lmtp::lmtp_tmle(outcome_type = "survival")`,
  `gfoRmula::gformula_survival()`, `survival::survfit` /
  `survival::coxph`, or a closed-form DGP).
- 🟡 — smoke-tested only (runs without error; point estimate / SE not
  pinned to a truth). Acceptable for combinations where no oracle exists,
  temporary during a multi-chunk feature rollout, or where the oracle is
  too expensive to run in CI.
- 🔴 — hard-rejected by a classed error with a regression test pinning the
  rejection.

<!-- Rows are added as features ship. No speculative ⚪ entries: the matrix
reflects **current** state, not planned scope. Planned scope lives in
`SURVIVAL_PACKAGE_HANDOFF.md` §10 (implementation chunks). -->
