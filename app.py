import pandas as pd
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import HTML
import numpy as np
import json
from multiprocessing import Pool

from Bio import SeqIO
from Bio import Align

import json
import numpy as np
from functools import reduce
from utils import load_patients, load_viruses, compute_virus_alignment_score

plt.ioff()

HTML("""\
<style>
.app-subtitle {
    font-size: 1.5em;
}

.app-subtitle a {
    color: #106ba3;
}

.app-subtitle a:hover {
    text-decoration: underline;
}

.app-sidebar p {
    margin-bottom: 1em;
    line-height: 1.7;
}

.app-sidebar a {
    color: #106ba3;
}

.app-sidebar a:hover {
    text-decoration: underline;
}
</style>
""")


class App:

    def __init__(self):
        with open('popular_virus_names.json', 'r') as f:
            popular_viruses_keys = json.load(f)

        self._virus_dropdown = self._create_indicator_dropdown(popular_viruses_keys, 0)
        self._plot_container = widgets.Output()
        _app_container = widgets.VBox([
            widgets.HBox([self._virus_dropdown]),
            self._plot_container,
        ], layout=widgets.Layout(align_items='center', flex='3 0 auto'))
        self.container = widgets.VBox([
            widgets.HTML(
                (
                    '<h1>Development indicators</h1>'
                    '<h1> </h1>'
                    '<h3>Select virus of interest</h3>'
                ),
                layout=widgets.Layout(margin='0 0 2em 0')
            ),
            widgets.HBox([
                _app_container,
            ])
        ], layout=widgets.Layout(flex='1 1 auto', margin='0 auto 0 auto', max_width='1024px'))

        # Loading and preprocessing data.
        n_viruses = 5
        n_entries = 10
        # n_viruses = 87
        # n_entries = 100

        viruses_full = load_viruses()
        patients_identity_full, patients_entries_full = load_patients()
        self.viruses = {k: viruses_full[k] for k in popular_viruses_keys[:n_viruses]}
        random_idx = np.arange(patients_entries_full.shape[0])
        np.random.shuffle(random_idx)
        random_idx = random_idx[:n_entries]
        patients_entries = patients_entries_full[random_idx]
        self.patients_identity = patients_identity_full[random_idx]
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -1
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0

        def comp_al_score(patients_vdna):
            patients_entries, virus_dna = patients_vdna
            return compute_virus_alignment_score(patients_entries, virus_dna, aligner)

        tuples_to_process = [(patients_entries, v_dna) for v_dna in self.viruses.values()]
        # with Pool(8) as p:
        #     alignements_scors_list = p.map(
        #         comp_al_score,
        #         tuples_to_process
        #     )
        alignements_scors_list = list(map(lambda x: comp_al_score(x), tuples_to_process))
        self.result = np.array(alignements_scors_list).T
        patients_results = dict()
        for patient_id, entry in zip(self.patients_identity, self.result):
            if patient_id not in patients_results:
                patients_results[patient_id] = []
            else:
                patients_results[patient_id] += [entry]
        for k in patients_results.keys():
            patients_results[k] = np.array(patients_results[k])
        patients_means_per_virus = dict()
        for k in patients_results.keys():
            patients_means_per_virus[k] = np.median(patients_results[k], axis=0)

    def _plot_histograms(self, virus_to_look):
        result_copy = self.result.copy()

        patients_results = dict()
        for patient_id, entry in zip(self.patients_identity, result_copy):
            if patient_id not in patients_results:
                patients_results[patient_id] = []
            else:
                patients_results[patient_id] += [entry]

        for k in patients_results.keys():
            patients_results[k] = np.array(patients_results[k])

        patients_means_per_virus = dict()
        for k in patients_results.keys():
            patients_means_per_virus[k] = np.median(patients_results[k], axis=0)

        quantiles_label = [0.25, 0.5, 0.75, 1]
        virus_quantiles = np.quantile(result_copy, q=quantiles_label, axis=0)

        patient_keys = list(patients_means_per_virus.keys())
        virus_idx = list(self.viruses.keys()).index(virus_to_look)
        patients_points = list()
        for k in patient_keys:
            if patients_means_per_virus[k] is not None and not np.isnan(patients_means_per_virus[k]).all():
                patients_points.append(patients_means_per_virus[k][virus_idx])
            else:
                patients_points.append(0)

        plt.figure(num=None, figsize=(15, 10), dpi=250, facecolor='w', edgecolor='k')
        min_max_x = (-0.05, len(patients_points) - 1 + 0.05)

        last_v = 0
        for name, v, c in zip(quantiles_label, virus_quantiles[:, virus_idx], ['green', 'yellow', 'orange', 'red']):
            plt.axhline(y=v, linestyle='--', label=f'{float(name):.2%} quartile', color=c)
            plt.fill_between(x=min_max_x, y1=last_v, y2=v, color=c, alpha=0.6)
            last_v = v
        plt.gca().tick_params(axis='both', which='major', labelsize=16)
        plt.scatter(x=patient_keys, y=patients_points, s=300, marker='x', c='black')
        for x, (name, y) in enumerate(zip(patient_keys, patients_points)):
            plt.annotate(name, (x + 0.05, y))
        plt.xlim(*min_max_x)
        plt.legend()

    def _create_indicator_dropdown(self, indicators, initial_index):
        dropdown = widgets.Dropdown(options=indicators, value=indicators[initial_index])
        dropdown.observe(self._on_change, names=['value'])
        return dropdown

    def _on_change(self, _):
        self._update_app()

    def _update_app(self):
        virus_indicator = self._virus_dropdown.value
        self._plot_container.clear_output(wait=True)
        with self._plot_container:
            self._plot_histograms(virus_indicator)
            plt.show()
