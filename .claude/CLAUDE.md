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
2. Group the existing open PRs into a stack on GitHub by listing them
   bottom → top (PR numbers, branches, or URLs):
   ```
   gh stack link <bottom-PR> <next-PR> ... <top-PR>
   ```
   `gh stack link` adopts the existing open PRs — it does not recreate branches.
3. Adopt the stack locally — `link` writes GitHub state only, so until then
   `gh stack view` reports "not part of a stack":
   ```
   gh stack checkout <stack-number>    # the number link printed
   ```
4. Inspect the stack with `gh stack view`, from a branch that is part of it.

Do NOT run `gh stack init` to adopt existing branches: in the installed
extension version it re-creates the listed branches from the current HEAD and
collapses the whole stack onto one commit. Use `gh stack link` instead.

PR branches follow the naming convention: `feat/acoustics-<slug>`

### Rebasing a stack

`gh stack rebase` cascades across the whole stack, so a conflict is resolved with
the stack commands, not the plain git ones:

```
# resolve the conflicted files, then
git add <resolved-files>
gh stack rebase --continue     # NOT git rebase --continue
gh stack rebase --abort        # to back out the whole cascade
```

`git rebase --continue` finishes only the branch that stopped and leaves the
remaining layers unrebased. Note that continuing writes a commit, so ask before
running it.

## Commit messages and PR descriptions

Describe only the delta relative to the parent branch. Never describe a change as
a fix to code introduced earlier in the same branch — that code does not exist
from the parent's point of view. Fold the correction into the feature
description instead.
