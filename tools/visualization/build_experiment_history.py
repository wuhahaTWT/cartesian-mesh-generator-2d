#!/usr/bin/env python3
"""Build a read-only history atlas. No solver runs, uploads, or deletion.

File dispositions are screening suggestions, not permission to delete. Historical
claims come from CURRENT_STATE_CN.md; file timestamps are not experiment dates.
"""
import collections
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'artifacts/experiment-history'
GROUPS = ['网格与质量', '桌面与输入输出', '外部 CFD 对照', '原生层流与精度',
          '非定常与续算', '标量与传热', '求解性能', 'SST 湍流', '范围与归档']


def topic(title):
    if any(x in title for x in ['SST', '平板', '近壁', '薄矩形', '高雷诺', '高低位', '空间黏度']): return 'SST 湍流'
    if any(x in title for x in ['热', '标量', '空间扩散']): return '标量与传热'
    if any(x in title for x in ['非定常', '物理时间', '时间步']): return '非定常与续算'
    if any(x in title for x in ['性能', '稀疏', 'IC(0)', '多重网格', '动量装配', '压力修正', '残差范数']): return '求解性能'
    if any(x in title for x in ['OpenFOAM', '已有 CFD']): return '外部 CFD 对照'
    if any(x in title for x in ['方腔', '流动', '层流', 'FVM', '压力', '应力']): return '原生层流与精度'
    if any(x in title for x in ['桌面', 'App', '图片', 'PNG', '输入', '导入']): return '桌面与输入输出'
    if any(x in title for x in ['核心', '质量', '密度', '缩放', '规模证据']): return '网格与质量'
    return '范围与归档'


def case_for(rel):
    parts = rel.split('/')
    if parts[0] == 'native-flow' and len(parts) > 2: return '/'.join(parts[:2])
    return parts[0] if len(parts) > 1 else '[根目录文件]'


def disposition(rel, size):
    p = Path(rel); ext = p.suffix.lower(); name = p.name
    if rel.startswith('archive-selection/'):
        return '保留记录', '归档索引和恢复记录；迁出 outputs 后保留'
    if '.app/' in rel or ext in ['.o', '.a', '.dylib', '.pak', '.pyc'] or name == '.DS_Store':
        return '可重建候选', '旧程序、构建产物或缓存；确认无唯一用途后清理'
    if ext in ['.png', '.jpg', '.jpeg', '.svg']:
        return '保留记录', '历史图像；精选展示图，其余可压缩归档'
    if ext in ['.py', '.sh', '.cpp', '.hpp', '.cjs', '.js', '.md', '.xy', '.dxf', '.geo']:
        return '保留记录', '输入或复现材料；先核对是否已有源码版本'
    if size <= 262144 and ext in ['.json', '.txt', '.log', '.csv', '.dat', '.jsonl']:
        return '保留记录', '小型记录候选；不自动把所有 JSON 当作结论'
    if name in ['controlDict', 'fvSchemes', 'fvSolution', 'transportProperties', 'turbulenceProperties', 'boundary', 'meshQualityDict']:
        return '保留记录', '工况或边界配置'
    if ext in ['.cm2d', '.vtk', '.vtp', '.checkpoint', '.csv', '.json'] or name in ['U', 'p', 'phi', 'k', 'omega', 'nut', 'points', 'faces', 'owner', 'neighbour', 'T']:
        return '结果待精选', '完整网格、场或预览；里程碑保留一份，其余在结论提取后清理'
    return '用途待核对', '仅凭名称不能判定；需确认唯一输入、结果或失败案例'


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'docs/CURRENT_STATE_CN.md'
    source_text = source.read_text()
    lines = source_text.splitlines()
    starts = [i for i, line in enumerate(lines) if line.startswith('## ')]
    sections = []
    for i, start in enumerate(starts):
        end = starts[i + 1] if i + 1 < len(starts) else len(lines)
        body = '\n'.join(lines[start+1:end]).strip()
        title = lines[start][3:]
        evidence = []
        for match in re.finditer(r'\[([^\]]+)\]\(([^)]+)\)', body):
            label, link = match.groups()
            if not link.startswith('../artifacts/'): continue
            path = (source.parent / link).resolve()
            if ROOT not in path.parents: continue
            evidence.append({'label': label, 'path': str(path.relative_to(ROOT)), 'exists': path.is_file(), 'image': path.suffix.lower() in ['.png', '.jpg']})
        sections.append({'id': i, 'title': title, 'topic': topic(title), 'line': start+1,
                         'body': body, 'evidence': evidence,
                         'abstract': next((p for p in body.split('\n\n') if p.strip()), '')})
    index_file = DEST / 'evidence-index.json'
    if index_file.is_file():
        extra_evidence = json.loads(index_file.read_text())
        for section in sections:
            known = {entry['path'] for entry in section['evidence']}
            for entry in extra_evidence.get(section['title'], []):
                path = (ROOT / entry['path']).resolve()
                if ROOT not in path.parents:
                    raise ValueError('Evidence path outside repository')
                if entry['path'] not in known:
                    section['evidence'].append(dict(entry, exists=path.is_file(), image=path.suffix.lower() in ['.png', '.jpg']))
                    known.add(entry['path'])
    records = []; symlinks = []; errors = []
    for base, ds, fs in os.walk(ROOT / 'outputs', followlinks=False, onerror=lambda e: errors.append(str(e))):
        ds.sort(); fs.sort()
        for d in list(ds):
            p = Path(base) / d
            if p.is_symlink(): symlinks.append(str(p.relative_to(ROOT))); ds.remove(d)
        for f in fs:
            p = Path(base) / f
            if p.is_symlink(): symlinks.append(str(p.relative_to(ROOT))); continue
            try: st = p.stat()
            except OSError as e: errors.append(str(e)); continue
            rel = str(p.relative_to(ROOT / 'outputs'))
            action, reason = disposition(rel, st.st_size)
            records.append({'path': 'outputs/' + rel, 'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns,
                            'case': case_for(rel), 'action': action, 'reason': reason})
    if errors: raise RuntimeError('\n'.join(errors))
    cases = {}
    for r in records:
        c = cases.setdefault(r['case'], {'name': r['case'], 'bytes': 0, 'files': 0, 'actions': collections.Counter(), 'largest': [], 'sections': []})
        c['bytes'] += r['bytes']; c['files'] += 1; c['actions'][r['action']] += r['bytes']
        c['largest'].append([r['bytes'], r['path'], r['action']])
    for c in cases.values():
        c['largest'] = sorted(c['largest'], reverse=True)[:8]
        needle = 'outputs/' + c['name']
        c['sections'] = [s['id'] for s in sections if re.search(re.escape(needle) + r'(?=[/`\s)。，；]|$)', s['body'])]
    rawlog = subprocess.check_output(['git', 'log', '--reverse', '--format=%h%x09%aI%x09%s'], cwd=ROOT, text=True)
    commits = [dict(zip(['sha', 'date', 'title'], l.split('\t', 2))) for l in rawlog.splitlines()]
    totals = collections.Counter()
    for r in records: totals[r['action']] += r['bytes']
    data = {'generated': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
            'commit': subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, text=True).strip(),
            'bytes': sum(r['bytes'] for r in records), 'files': len(records), 'sections': sections,
            'cases': sorted(cases.values(), key=lambda c: -c['bytes']), 'commits': commits,
            'actions': dict(totals), 'topics': GROUPS, 'symlinks': symlinks,
            'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest()}
    inventory = json.dumps({'scope': 'Proposed screening only; no deletion authorized by this file; no cloud verification performed', 'generated': data['generated'], 'files': records, 'symlinks_not_followed': symlinks}, ensure_ascii=False, separators=(',', ':')).encode()
    (DEST / 'file-inventory.json.gz').write_bytes(gzip.compress(inventory, mtime=0))
    (DEST / 'history-data.json').write_text(json.dumps(data, ensure_ascii=False, separators=(',', ':')))
    for field, filename in [('cleanup', 'cleanup-result.json'), ('packageCleanup', 'package-cleanup-result.json'), ('seafileReview', 'seafile-review.json')]:
        if (DEST / filename).is_file():
            data[field] = json.loads((DEST / filename).read_text())
    (DEST / 'history-data.json').write_text(json.dumps(data, ensure_ascii=False, separators=(',', ':')))
    template = Path(__file__).with_name('experiment_history_template.html').read_text()
    payload = json.dumps(data, ensure_ascii=False, separators=(',', ':')).replace('<', '\\u003c')
    (DEST / 'index.html').write_text(template.replace('__HISTORY_DATA__', payload))
    print(json.dumps({'files': len(records), 'bytes': data['bytes'], 'cases': len(cases), 'document_sections': len(sections), 'explicitly_linked_cases': sum(bool(c['sections']) for c in cases.values()), 'commits': len(commits), 'actions': totals, 'html_bytes': (DEST/'index.html').stat().st_size}, ensure_ascii=False))


if __name__ == '__main__': main()
