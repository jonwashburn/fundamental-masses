import Lake
open Lake DSL

package «ParticleMasses» where
  leanOptions := #[
    ⟨`autoImplicit, false⟩
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git"

@[default_target]
lean_lib «IndisputableMonolith» where
  -- add library configuration options here
