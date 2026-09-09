"""Independent reader regression: real wall partition with a stopped column."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.verification.check_layer_resolution import CheckError, measure


class StoppedColumnTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.vtk, self.cm2d, self.report = (root / name for name in ('source.vtk', 'final.cm2d', 'resolution.json'))
        points = [(x, y) for y in range(3) for x in range(3)]
        cells = [(0,1,4,3), (1,2,5,4), (3,4,7,6), (4,5,8,7)]
        edge_map, edges, cell_edges = {}, [], []
        for cid, vertices in enumerate(cells):
            ids = []
            for a, b in zip(vertices, vertices[1:]+vertices[:1]):
                key = tuple(sorted((a,b)))
                if key in edge_map:
                    eid = edge_map[key]
                    edges[eid][3] = cid
                    edges[eid][4] = 0
                else:
                    eid = len(edges); edge_map[key] = eid
                    patch = 1 if points[a][1] == points[b][1] == 0 else 2
                    edges.append([a,b,cid,-1,patch])
                ids.append(eid)
            cell_edges.append(ids)
        text = ['CM2D 1', f'VERTICES {len(points)}']
        text += [f'{i} {x} {y}' for i,(x,y) in enumerate(points)]
        text += [f'EDGES {len(edges)}']
        text += [' '.join(map(str, [i]+edge)) for i,edge in enumerate(edges)]
        text += [f'CELLS {len(cells)}']
        text += [' '.join(map(str,[i,i,i,1,4,*vertices,4,*cell_edges[i]])) for i,vertices in enumerate(cells)]
        text += ['AUDIT 0 0 0 0 0 0 0', 'END']
        self.cm2d.write_text('\n'.join(text)+'\n')
        text = ['# vtk DataFile Version 3.0', 'four-cell stopped-layer fixture', 'ASCII',
                'DATASET UNSTRUCTURED_GRID', f'POINTS {len(points)} double']
        text += [f'{x} {y} 0' for x,y in points]
        text += ['CELLS 4 20'] + ['4 '+' '.join(map(str,c)) for c in cells]
        text += ['CELL_TYPES 4', '9 9 9 9', 'CELL_DATA 4',
                 'SCALARS hybrid_kind int 1', 'LOOKUP_TABLE default', '0 2 2 2',
                 'SCALARS layer_index int 1', 'LOOKUP_TABLE default', '0 -1 -1 -1']
        self.vtk.write_text('\n'.join(text)+'\n')
        columns = [dict(strip_id=0,wall_segment=i,requested_layers=1,retained_layers=1-i,
                        wall_length_over_reference=1,wall_segment_endpoints=[[i,0],[i+1,0]],
                        solver_cell_ids=[0] if i==0 else [],
                        first_layer_normal_height_min_over_reference=1 if i==0 else None,
                        first_layer_normal_height_max_over_reference=1 if i==0 else None) for i in range(2)]
        self.data = dict(reference_length=1,requested=dict(first_layer_h_over_reference=1),
            boundary_layers=dict(status='evaluated',columns=columns,requested_cells=2,
                constructed_cells=1,retained_cells=1,source_mismatches=0,continuity_mismatches=0,
                first_layer_wall_length_fraction=.5,full_requested_layers_wall_length_fraction=.5,
                first_layer_height_exceedance_wall_length_fraction=0))

    def read(self, data=None):
        self.report.write_text(json.dumps(self.data if data is None else data))
        return measure(self.vtk,self.cm2d,self.report)

    def test_stopped_column_has_no_layer_polygon_but_real_wall(self):
        result = self.read()
        self.assertTrue(result['valid'])
        self.assertEqual(result['source_wall_length'],2)
        self.assertEqual(result['first_layer_wall_length_fraction'],.5)
        self.assertIsNone(result['columns'][1]['first_layer_height_max'])

    def test_reject_false_wall_coverage_and_height_claims(self):
        for defect in ('duplicate', 'missing', 'domain_side', 'gap', 'height', 'coverage', 'length'):
            with self.subTest(defect=defect):
                data = copy.deepcopy(self.data); boundary=data['boundary_layers']; columns=boundary['columns']
                if defect=='duplicate': columns.append(copy.deepcopy(columns[1]))
                elif defect=='missing': columns.pop()
                elif defect=='domain_side': columns[1]['wall_segment_endpoints']=[[0,0],[0,1]]
                elif defect=='gap': columns[1]['wall_segment_endpoints']=[[1.1,0],[2,0]]; columns[1]['wall_length_over_reference']=.9
                elif defect=='height': columns[1]['first_layer_normal_height_max_over_reference']=0
                elif defect=='coverage': boundary['first_layer_wall_length_fraction']=1
                elif defect=='length': columns[1]['wall_length_over_reference']=2
                with self.assertRaises((CheckError,ValueError,TypeError)):
                    self.read(data)

    def test_reject_layer_polygon_source_mismatch(self):
        text=self.vtk.read_text().replace('0 2 2 2','2 0 2 2')
        self.vtk.write_text(text)
        with self.assertRaises(CheckError): self.read()


if __name__=='__main__': unittest.main()
