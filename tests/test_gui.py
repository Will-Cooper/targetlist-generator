import pytest

pytest.importorskip('PySide6')
pytestmark = pytest.mark.gui


@pytest.fixture
def window(qtbot):
    from targetlist_generator.gui import MainWindow
    win = MainWindow()
    qtbot.addWidget(win)
    yield win
    win.samp.stop()


def test_gui_plan_prompts_for_folder(qtbot, window, target_csv, tmp_path, monkeypatch):
    from PySide6.QtCore import QDate
    from PySide6.QtWidgets import QFileDialog
    asked = []
    monkeypatch.setattr(QFileDialog, 'getExistingDirectory', lambda *a: asked.append(a[1]) or str(tmp_path))
    window.load_targets(str(target_csv))
    assert window.mag_col.currentText() == 'tmassk'
    assert window.site.currentText().startswith('NASA Infrared Telescope Facility')
    assert window.lat.value() == pytest.approx(19.826, abs=1e-3)
    window.date.setDate(QDate(2026, 11, 30))
    window.format_boxes['keck'].setChecked(True)
    window.format_boxes['csv'].setChecked(False)
    with qtbot.waitSignal(window.finished, timeout=120_000):
        window.run_button.click()
    assert asked == ['Where should the lists and charts be saved?']
    headers = [window.table.horizontalHeaderItem(j).text() for j in range(window.table.columnCount())]
    assert {'status', 'target_start', 'standard', 'std_sep'} <= set(headers)
    assert (tmp_path / 'targets_keck.txt').exists() and not (tmp_path / 'targets_list.csv').exists()
    assert window.plots.count() >= 1 and not window.plot.pixmap().isNull()


def test_bundled_qt_plugins(monkeypatch):
    """A conda qt6.conf must not send PySide6 to another Qt's plugins; explicit user paths win."""
    from pathlib import Path

    import PySide6

    from targetlist_generator.gui import use_bundled_qt_plugins
    monkeypatch.delenv('QT_PLUGIN_PATH', raising=False)
    monkeypatch.setenv('QT_QPA_PLATFORM_PLUGIN_PATH', '/my/plugins')
    use_bundled_qt_plugins()
    import os
    assert os.environ['QT_PLUGIN_PATH'] == str(Path(PySide6.__file__).parent / 'Qt' / 'plugins')
    assert os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] == '/my/plugins'


def test_gui_custom_coordinates(window):
    from targetlist_generator.gui import CUSTOM
    window.lat.setValue(-24.6)
    assert window.site.currentText() == CUSTOM and window.site.currentData() is None


def test_gui_needs_target_list(window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    shown = []
    monkeypatch.setattr(QMessageBox, 'warning', lambda *args: shown.append(args[-1]))
    window.run_button.click()
    assert shown and shown[0].startswith('Choose a target table first')


def test_gui_receives_samp_table(qtbot, window, target_csv, tmp_path):
    """A second SAMP client (standing in for TOPCAT) sends a VOTable; the window loads it."""
    from astropy.table import Table
    from pyvo.samp import SAMPHubServer, SAMPIntegratedClient
    hub = SAMPHubServer(web_profile=False, lockfile=str(tmp_path / 'samp.lock'))
    hub.start()
    import os
    os.environ['SAMP_HUB'] = f'std-lockurl:file://{tmp_path / "samp.lock"}'
    try:
        window.listen.setChecked(True)
        assert window.samp.client is not None and window.samp.hub is None  # joined the running hub
        vot = tmp_path / 'sent.vot'
        Table.read(target_csv).write(vot, format='votable')
        topcat = SAMPIntegratedClient(name='fake topcat')
        topcat.connect()
        try:
            with qtbot.waitSignal(window.samp.table, timeout=20_000):
                topcat.notify_all({'samp.mtype': 'table.load.votable',
                                   'samp.params': {'url': vot.as_uri(), 'name': 'sent'}})
            assert window.target_list.text().endswith('.vot') and window.mag_col.currentText() == 'tmassk'
        finally:
            topcat.disconnect()
    finally:
        window.listen.setChecked(False)
        os.environ.pop('SAMP_HUB', None)
        hub.stop()
