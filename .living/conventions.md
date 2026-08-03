# Repo-Specific Conventions

Overrides to mycelium defaults or convention pack conventions.

## Development Conventions

### Nextflow Schema Validation (nextflow_schema.json)

**Rule**: When making changes to any file under `nextflow/module/`, `nextflow/submodule/`, `nextflow/workflow/`, or `nextflow/conf/`, validate that `nextflow/nextflow_schema.json` remains up to date before committing.

**Why**: The schema defines parameter validation and documentation for the Nextflow pipeline. Failing to update it when adding/changing parameters leads to: (1) silent runtime failures if parameters are omitted, (2) incorrect validation of user inputs, (3) broken CLI help/documentation generation, (4) downstream bugs that are hard to trace back to schema drift.

**How to apply**:
1. After editing any Nextflow configuration, module, submodule, or workflow file, check `nextflow/nextflow_schema.json` for completeness.
2. If you added or renamed a parameter in `nextflow/conf/*.config`, ensure it's documented in the schema with type, description, and default (if applicable).
3. If you added a new input channel or parameter to a module/submodule, verify it doesn't need schema documentation (schema typically covers pipeline-level params, not module internals, but confirm this matches this repo's pattern).
4. Run `nextflow run ... -profile test --help` to confirm the schema is valid and help text renders correctly.
5. A schema validation error will surface as a runtime error; treat it as a blocker to merge.

**Related**: `nextflow/nextflow_schema.json`, [[project_nfschema_validation]].
