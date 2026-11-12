import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

from .evidence import build_evidence
from .io import load_records
from .summary import summarize


def save_json(path, result):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix='.analysis-', dir=destination.parent)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, sort_keys=True, allow_nan=False, indent=2)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def save_text(path, text):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix='.analysis-', dir=destination.parent)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as handle:
            handle.write(text)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description='本地观测数据分析')
    commands = parser.add_subparsers(dest='command', required=True)
    summary = commands.add_parser('summary', help='按业务维度核对观测口径')
    summary.add_argument('--input', type=Path, default=Path('data/observations.csv'))
    summary.add_argument('--output', type=Path)
    evidence = commands.add_parser('evidence', help='生成可复算的水色-浓度关联证据包')
    evidence.add_argument('--input', type=Path, default=Path('data/observations.csv'))
    evidence.add_argument('--cutoff', type=int, required=True, help='包含 day<=cutoff 的观测行')
    evidence.add_argument('--seed', type=int, required=True, help='站点聚类自助法固定种子')
    evidence.add_argument('--iterations', type=int, default=1000, help='站点聚类自助次数')
    evidence.add_argument('--output', type=Path, help='机器可读 JSON 输出路径')
    evidence.add_argument('--narrative', type=Path, help='中文说明输出路径（.txt/.md）')
    args = parser.parse_args()
    try:
        if args.command == 'summary':
            result = summarize(load_records(args.input))
        else:
            result = build_evidence(load_records(args.input), args.cutoff,
                                    args.seed, args.iterations)
            if args.narrative:
                save_text(args.narrative, result['narrative_zh'])
        if args.output:
            save_json(args.output, result)
        print(json.dumps(result, ensure_ascii=False, sort_keys=True, allow_nan=False))
    except (ValueError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
