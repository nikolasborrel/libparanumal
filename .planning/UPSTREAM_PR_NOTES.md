# Upstream PR Notes

## Excluding workflow files from upstream diffs

The `.planning/` and `.claude/` directories are GSD workflow artifacts and should not appear in PRs submitted to `paranumal/libParanumal`.

Before submitting Phase 1 PR to upstream, run:

```bash
git rm -r --cached .planning .claude
echo ".planning/" >> .gitignore
echo ".claude/" >> .gitignore
git commit -m "chore: exclude workflow planning docs from upstream"
```

This removes them from git tracking so they vanish from the upstream PR diff, while:
- Files still exist on disk locally
- Files are still pushed to your fork (`nikolasborrel/libparanumal`)
- GSD continues to work normally
- No future `.planning/` changes will appear in phase PR diffs

## PR strategy

- Submit PRs one-by-one (or as a stack) to `paranumal/libParanumal`
- Reach out to maintainers via a GitHub issue first to get buy-in before opening PRs
- Phase 1 (tests) goes first — use it as the conversation starter
- Each phase targets the previous phase branch (stacked), except Phase 1 which targets `main`

## Remotes

- `origin` → `nikolasborrel/libparanumal` (your fork, write access)
- `upstream` → `paranumal/libParanumal` (read only)

## gh stack notes

- After any out-of-band force push, run `gh stack rebase` before `gh stack push`
