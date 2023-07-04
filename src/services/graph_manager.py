import os
import matplotlib.pyplot as plt
import numpy as np

from config import LOG_FILE_DIR


class GraphManagement:

    def __init__(self):
        self.graph_dir = os.path.join(LOG_FILE_DIR, "img")
        if not os.path.exists(self.graph_dir):
            os.mkdir(self.graph_dir)

    def normalize_values(self, graph_values: dict):
        return {key: value / sum(graph_values.values()) for key, value in graph_values.items()}

    def construct_bar_chart_from_dict(self,
                                      graph_values: dict,
                                      filename: str,
                                      title: str,
                                      x_label: str,
                                      y_label: str):
        normalized_graph_values = self.normalize_values(graph_values)
        values = sorted(normalized_graph_values.items())
        x_values, y_values = zip(*values)

        plt.bar(x_values, height=y_values, tick_label=x_values)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title, wrap=True)
        plt.savefig(os.path.join(self.graph_dir, filename + ".png"))
        plt.cla()
        plt.clf()
        plt.close()


default_graph_manager = GraphManagement()
