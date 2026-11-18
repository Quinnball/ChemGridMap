"""Verify that text-only edits left Python token streams unchanged."""
from __future__ import annotations

import argparse
import io
import tokenize
from pathlib import Path


REMOVED_TOKENS = {
    tokenize.COMMENT,
    tokenize.STRING,
    tokenize.FSTRING_START,
    tokenize.FSTRING_MIDDLE,
    tokenize.FSTRING_END,
}


def code_tokens(path: Path) -> list[tuple[int, str]]:
    text = path.read_text(encoding="utf-8")
    return [
        (token.type, token.string)
        for token in tokenize.generate_tokens(io.StringIO(text).readline)
        if token.type not in REMOVED_TOKENS
        and token.type not in {tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT, tokenize.ENDMARKER}
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    failed = []
    for name in args.files:
        before = args.before / name
        after = args.after / name
        if code_tokens(before) != code_tokens(after):
            failed.append(name)
    for name in failed:
        print(f"TOKEN_DIFF {name}")
    print(f"checked={len(args.files)} failed={len(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
