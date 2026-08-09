# libparanumal — Claude Code project instructions

## Stacked PRs

Create branches and PRs with plain git + `gh pr create`, then group them into a
stack with `gh stack link`. Do NOT use `gh stack submit`.

1. Branch and push each PR, chaining each `--base` on the branch below:
   ```
   git checkout -b feat/<slug> <base-branch>
   # implement changes, commit
   git push -u origin feat/<slug>
   gh pr create --base <base-branch> --head feat/<slug> ...
   ```
2. Group the existing open PRs into a gh-stack (local + GitHub tracking) by
   listing them bottom → top (PR numbers, branches, or URLs):
   ```
   gh stack link <bottom-PR> <next-PR> ... <top-PR>
   ```
   `gh stack link` adopts the existing open PRs — it does not recreate branches.
3. Inspect the stack with `gh stack view`.

Do NOT run `gh stack init` to adopt existing branches: in the installed
extension version it re-creates the listed branches from the current HEAD and
collapses the whole stack onto one commit. Use `gh stack link` instead.

PR branches follow the naming convention: `feat/acoustics-<slug>`
