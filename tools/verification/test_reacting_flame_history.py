#!/usr/bin/env python3
"""Analytic kinematics checks; no combustion-accuracy thresholds or solver runs."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from analyze_reacting_flame_history import fit_motion, marker_positions, sample_history, translated_shape


class FlameHistoryKinematics(unittest.TestCase):
    def test_affine_isotherm_on_nonuniform_centres(self):
        x = np.array([0., .1, .25, .7, .85, 1.])
        speed, level = -.07, 800.
        for time in [0., .23, .81]:
            # This translated affine profile has an exact isotherm at .4+c*t.
            temperature = 400 + 1000 * (x - speed * time)
            expected = .4 + speed * time
            result = marker_positions(x, temperature, level)
            self.assertEqual(len(result['crossings_m']), 1)
            self.assertEqual(result['plateaus_m'], [])
            # Position allowance: 64 eps times the 1 m coordinate scale.
            self.assertLessEqual(abs(result['unique_position_m'] - expected), 64 * np.finfo(float).eps * np.max(abs(x)))

    def test_sample_at_exact_crossing_is_counted_once(self):
        result = marker_positions([0., .2, .9], [400., 1000., 2000.], 1000.)
        self.assertEqual(result, {'crossings_m': [.2], 'plateaus_m': [], 'unique_position_m': .2})

    def test_two_crossings_cannot_be_selected_silently(self):
        result = marker_positions([0., .2, 1.], [400., 1200., 400.], 800.)
        self.assertEqual(len(result['crossings_m']), 2)
        self.assertLessEqual(np.max(abs(np.array(result['crossings_m']) - [.1, .6])), 64 * np.finfo(float).eps)
        self.assertIsNone(result['unique_position_m'])
        self.assertFalse(fit_motion([0., 1., 2.], [None, None, None])['available'])

    def test_plateau_and_absent_marker_are_explicit(self):
        plateau = marker_positions([0., .1, .3, .7, 1.], [400., 800., 800., 800., 1000.], 800.)
        self.assertEqual(plateau['plateaus_m'], [[.1, .7]])
        self.assertIsNone(plateau['unique_position_m'])
        absent = marker_positions([0., 1.], [400., 700.], 800.)
        self.assertEqual(absent['crossings_m'], [])
        self.assertIsNone(absent['unique_position_m'])

    def test_unequal_microsecond_samples_preserve_velocity_units(self):
        times = np.array([2e-6, 9e-6, 21e-6, 80e-6])
        speed = -.07
        positions = .014 + speed * times
        result = fit_motion(times, positions)
        self.assertTrue(result['available'])
        # Velocity allowance propagates coordinate rounding over actual seconds.
        allowance = 64 * np.finfo(float).eps * np.max(abs(positions)) / (times[-1] - times[0])
        for key in ['lab_velocity_m_per_s', 'endpoint_secant_velocity_m_per_s']:
            self.assertLessEqual(abs(result[key] - speed), allowance)
        self.assertLessEqual(result['maximum_linear_fit_position_deviation_m'], 64 * np.finfo(float).eps * np.max(abs(positions)))

    def test_nonuniform_speed_is_not_reported_as_steady(self):
        times = np.array([0., .1, .4, 1.])
        result = fit_motion(times, .014 + .1 * times**2)
        self.assertTrue(result['available'])
        self.assertGreater(result['maximum_linear_fit_position_deviation_m'], 0.)
        self.assertGreater(result['interval_lab_velocity_max_m_per_s'], result['interval_lab_velocity_min_m_per_s'])
        self.assertNotIn('steady', result)

    def test_insufficient_and_unordered_times(self):
        self.assertFalse(fit_motion([0., 1.], [.1, .2])['available'])
        with self.assertRaises(AssertionError):
            fit_motion([0., 1., 1.], [.1, .2, .3])

    def test_exact_translation_uses_only_supported_overlap(self):
        x = np.array([0., .05, .2, .4, .8, 1.])
        initial = 400 + 1800 * x
        displacement = .12
        current = 400 + 1800 * (x - displacement)
        widths = np.array([.02, .08, .16, .24, .3, .2])
        result = translated_shape(x, initial, current, displacement, widths)
        self.assertEqual(result['excluded_columns'], 2)
        self.assertEqual(result['evaluated_columns'], 4)
        self.assertLess(result['evaluated_width_fraction'], 1.)
        # Temperature allowance: 64 eps times the largest constructed temperature.
        self.assertLessEqual(result['maximum_K'], 64 * np.finfo(float).eps * np.max(abs(initial)))
        deformed = translated_shape(x, initial, current + 7., displacement, widths)
        self.assertLessEqual(abs(deformed['maximum_K'] - 7.), 64 * np.finfo(float).eps * np.max(abs(initial)))
        self.assertLessEqual(abs(deformed['width_weighted_RMS_K'] - 7.), 64 * np.finfo(float).eps * np.max(abs(initial)))
        self.assertIsNone(translated_shape(x, initial, current, None, widths))
        self.assertIsNone(translated_shape(x, initial, current, 2., widths))

    def test_stream_requires_both_endpoints_count_and_monotone_time(self):
        # Minimal records exercise JSONL completeness only, not CFD validation.
        records = [{'state': {'time': t, 'steps': i}} for i, t in enumerate([0., .1, .3])]
        initial, final = records[0]['state'], records[-1]['state']
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'samples.jsonl'
            path.write_text(''.join(json.dumps(r) + '\n' for r in records))
            self.assertEqual(list(sample_history(path, initial, final, 3)), records)
            for invalid in [records[:-1], records + [{'state': {'time': .4, 'steps': 3}}],
                            [records[0], records[1], records[1]], records[1:]]:
                path.write_text(''.join(json.dumps(r) + '\n' for r in invalid))
                with self.assertRaises(AssertionError):
                    list(sample_history(path, initial, final, 3))


if __name__ == '__main__':
    unittest.main()
