# libparanumal — Claude Code project instructions

## Stacked PRs

Always use `gh stack` when creating stacked PRs for this project.
Do NOT use `gh stack submit` — no beta access. Create PRs manually with `gh pr create`.

Workflow:
1. From `main`: create the first feature branch and initialize the stack
   ```
   git checkout -b feat/<slug>
   gh stack init --adopt feat/<slug>
   # implement changes, commit
   gh stack push   # or: git push -u origin feat/<slug>
   gh pr create --base main --head feat/<slug> ...
   ```
2. For each subsequent PR in the stack:
   ```
   gh stack add feat/<next-slug>
   # implement changes, commit
   gh stack push
   gh pr create --base feat/<prev-slug> --head feat/<next-slug> ...
   ```
3. When rebasing after upstream changes:
   ```
   gh stack rebase
   ```

PR branches follow the naming convention: `feat/acoustics-<slug>`
