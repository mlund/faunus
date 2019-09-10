# Description

Please edit this text to include a summary of the change and which issue is fixed.

## Checklist

- [ ] `make test` passes with no errors
- [ ] the source code is well documented
- [ ] new functionality includes unittests
- [ ] the user-manual has been updated to reflect the changes
- [ ] I have installed the provided commit hook to auto-format commits (requires `clang-format`):

  ``` bash
  ./scripts/git-pre-commit-format install
  ```

- [ ] code naming scheme follows the convention:

  Style      | Elements
  ---------- | -------------------
  PascalCase | classes, namespaces
  camelCase  | functions
  snake_case | variables
