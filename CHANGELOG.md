# v0.2.0
- Update to Agents v6
- General API improvements and re-export of Agents.jl API functions
- Reexport `StaticVector`

## BREAKING
- `pos` and `vel` fields of microbe types now use `SVector` rather than `NTuple`
- `id`, `pos`, `vel` and `speed` are not given default values when constructed outside of a model
- `rand_vel` replaced by `random_velocity`
