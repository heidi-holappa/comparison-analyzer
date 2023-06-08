from unittest import TestCase
from unittest.mock import MagicMock, patch

from services.db_initializer import init_databases
from services import db_initializer


class TestCompana(TestCase):

    @patch('services.db_initializer.init_databases')
    def test_init_databases_is_called(self, mock_init_databases):
        mock_init_databases.return_value = MagicMock(), MagicMock()
        db_initializer.init_databases(MagicMock())
        mock_init_databases.assert_called_once()
