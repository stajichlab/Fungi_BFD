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

### CSV and TSV Writers Must Use Unix Line Endings

**Rule**: All Python scripts that write CSV or TSV files must explicitly specify `lineterminator='\n'` when creating `csv.writer()` or `csv.DictWriter()` objects.

**Why**: Line ending inconsistencies cause merge conflicts, diffs that show every line as changed even when content is identical, and platform-specific behavior that breaks portability. Unix line endings (`\n`) are the repository standard and required for clean git history.

**How to apply**:
1. When using `csv.writer()`, pass `lineterminator='\n'`:
   ```python
   writer = csv.writer(f, lineterminator='\n')
   ```
2. When using `csv.DictWriter()`, pass `lineterminator='\n'`:
   ```python
   writer = csv.DictWriter(f, fieldnames=header, lineterminator='\n')
   ```
3. For pandas, use `line_terminator='\n'`:
   ```python
   df.to_csv(path, line_terminator='\n', index=False)
   ```
4. Do NOT rely on OS defaults or Python's cross-platform `newline=''` parameter alone — explicitly set the terminator.

**Related**: Data curation scripts, rescue scripts in `scripts/`.
