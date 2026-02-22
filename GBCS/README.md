# GBCS Merge Planning

This folder stores merge planning notes and scripts for building a single cohort-wide main dataset from:

- `/Users/linxu/Documents/GBCS/fu_may2014.dta` (anchor cohort; baseline + first follow-up content)
- `/Users/linxu/Documents/GBCS/gbcs3+端粒.dta` (second repeated measurement / wave 2)
- `/Users/linxu/Documents/GBCS/gbcs4+端粒.dta` (third repeated measurement / wave 3)

Current scope:

- Planning only (no merge execution in this step)
- Duplicate-resolution rule defined for later implementation
- Git commits/pushes include only code and markdown documentation

Data handling guardrails:

- Do not commit `.dta` source files
- Do not commit merged outputs (`.csv`, `.rds`, `.RData`, `.sav`)
- Use `obje_id` as the canonical person-level merge key
- Treat `obje_id_2f` / `obje_id_3f` as visit-specific IDs for audit/QC
