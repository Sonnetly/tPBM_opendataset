"""Audit actual measured MERFISH panels using official HTML gene lists."""
from pathlib import Path
from html.parser import HTMLParser
import hashlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


class GeneTable(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.row, self.cell = [], None, None

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.row = []
        if tag in {'th', 'td'}:
            self.cell = []

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag in {'th', 'td'} and self.cell is not None:
            self.row.append(''.join(self.cell).strip())
            self.cell = None
        if tag == 'tr' and self.row:
            self.rows.append(self.row)
            self.row = None


def read_table(path):
    parser = GeneTable()
    parser.feed(path.read_text(encoding='utf-8'))
    return pd.DataFrame(parser.rows[1:], columns=parser.rows[0])


def main():
    output = ROOT / 'results' / 'coverage'
    output.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(ROOT / 'results/pbm/panel22.csv')
    exploration = pd.read_csv(ROOT / 'results/pbm/exploration_genes62.csv')
    inputs = {'B6J': ('B6J_gene_list.html', 500),
              'B6N_Zhuang1': ('B6N_Zhuang1_gene_list.html', 1122)}
    results, sources = {}, {}
    for name, (file, count) in inputs.items():
        path = ROOT / 'input/coverage' / file
        table = read_table(path)
        assert len(table) == count, (name, len(table), count)
        assert not table.gene_symbol.duplicated().any(), name
        table.to_csv(output / f'{name}_measured_gene_list.csv', index=False, encoding='utf-8-sig')
        genes = set(table.gene_symbol)
        for data in (panel, exploration):
            data[name + '_measurement_status'] = data.gene.map(
                lambda g: 'measured_panel' if g in genes else 'unmeasured_not_zero')
        results[name] = {'panel_size': count,
                        'panel22_measured': sorted(set(panel.gene) & genes),
                        'exploration62_measured': sorted(set(exploration.gene) & genes)}
        sources[name] = {'relative_file': str(path.relative_to(ROOT)),
                         'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    assert results['B6J']['panel22_measured'] == ['Penk']
    assert results['B6N_Zhuang1']['panel22_measured'] == ['Adora1', 'Adora2a', 'Penk', 'Trpc5']
    panel.to_csv(output / 'panel22_measurement_coverage.csv', index=False, encoding='utf-8-sig')
    exploration.to_csv(output / 'exploration62_measurement_coverage.csv', index=False, encoding='utf-8-sig')
    summary = {'status': 'completed_measured_panel_audit', 'regional_expression_comparison': 'not_run',
               'results': results, 'sources': sources,
               'official_source_repository': 'https://github.com/AllenInstitute/abc_atlas_access',
               'source_commit': 'b16ee223d3466a57b30066973413ccc948c2e983',
               'source_paths': ['docs/gene_lists/MERFISH-C57BL6J-638850/gene_list.html',
                                'docs/gene_lists/Zhuang-ABCA-1/gene_list.html'],
               'interpretation': 'Panel absence means unmeasured; no zeros are imputed.'}
    (output / 'validation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
