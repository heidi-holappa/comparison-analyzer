import os
from unittest import TestCase

from services.graph_manager import default_graph_manager as graph_manager
from config import LOG_FILE_DIR


class TestGraphManager(TestCase):

    def test_graph_manager_creates_an_output_file(self):
        graphpath = os.path.join(LOG_FILE_DIR, "img")
        filepath = os.path.join(graphpath, "test_graph.png")
        if os.path.exists(filepath):
            os.remove(filepath)
        graph_values = {'a': 1, 'b': 2, 'c': 3}
        title = "test_graph"
        x_label = "x"
        y_label = "y"
        graph_manager.construct_bar_chart_from_dict(
            graph_values, "test_graph", title, x_label, y_label)
        assert os.path.exists(filepath)
        if os.path.exists(filepath):
            os.remove(filepath)
