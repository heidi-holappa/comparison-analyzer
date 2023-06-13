import os
import matplotlib.pyplot as plt
import numpy as np

from config import LOG_FILE_DIR


class GraphManagement:

    def __init__(self):
        self.graph_dir = os.path.join(LOG_FILE_DIR, "img")
        if not os.path.exists(self.graph_dir):
            os.mkdir(self.graph_dir)

    def construct_bar_chart_from_dict(self, graph_values: dict, title: str, x_label: str, y_label: str):
        values = sorted(graph_values.items())
        x_values, y_values = zip(*values)

        plt.bar(x_values, height=y_values, tick_label=x_values)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.savefig(os.path.join(self.graph_dir, title + ".png"))
        plt.cla()
        plt.clf()
        plt.close()
        # keys = list(graph_values.keys())
        # values = list(graph_values.values())

        # plt.bar(range(len(graph_values)), values, tick_label=keys)
        # plt.xlabel(x_label)
        # plt.ylabel(y_label)
        # plt.title(title)
        # plt.savefig(os.path.join(self.graph_dir, title + ".png"))


default_graph_manager = GraphManagement()
