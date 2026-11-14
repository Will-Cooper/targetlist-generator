import os
from pathlib import Path

import pandas as pd
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
    assert window.target_cols.choice()['mag'] == 'tmassk' and window.target_cols.choice()['shortname'] == 'name'
    assert window.table_status.text().startswith('25 targets') and window.run_button.isEnabled()
    assert window.site.currentText().startswith('NASA Infrared Telescope Facility')
    assert window.lat.value() == pytest.approx(19.826, abs=1e-3)
    window.date.setDate(QDate(2026, 11, 30))
    with qtbot.waitSignal(window.finished, timeout=120_000):
        window.run_button.click()
    assert asked == ['Where should the lists and charts be saved?']
    headers = [window.table.horizontalHeaderItem(j).text() for j in range(window.table.columnCount())]
    assert {'step', 'kind', 'start_utc', 'end_utc', 'airmass', 'separation_deg'} <= set(headers)
    kinds = [window.table.item(i, headers.index('kind')).text() for i in range(window.table.rowCount())]
    assert kinds[:2] == ['target', 'standard'] and kinds.count('target') == kinds.count('standard')
    assert (tmp_path / 'targets_irtf.txt').exists() and not (tmp_path / 'targets_list.csv').exists()
    assert len(window.items['visibility']) == window.plots.count() >= 1 and window.plot_view.count() == 1
    assert window.tabs.tabText(1) == f'Visibility ({window.plots.count()})'
    assert window.status.text().startswith('Saved') and window.progress.value() == 1


def test_default_formats_follow_observatory(window):
    from targetlist_generator.gui import CUSTOM
    ticked = lambda: sorted(k for k, box in window.format_boxes.items() if box.isChecked())  # noqa: E731
    assert ticked() == ['irtf']
    window.site.setCurrentIndex(next(i for i in range(window.site.count())
                                     if window.site.itemText(i).startswith('W. M. Keck')))
    assert ticked() == ['keck']
    window.site.setCurrentIndex(0)
    assert window.site.currentText() == CUSTOM
    assert ticked() == ['csv']
    window.from_samp = True
    window.default_formats()
    assert ticked() == ['votable']
    window.format_boxes['eso'].click()  # a user choice now stands
    window.site.setCurrentIndex(1)
    assert ticked() == ['eso', 'votable']


def test_columns_guessed_and_chosen(window, tmp_path):
    """Unrecognised names leave a field to choose (planning blocked); choosing it unblocks."""
    path = tmp_path / 'odd.csv'
    pd.DataFrame({'Designation': ['a'], 'alpha_deg': [10.0], 'DEJ2000': [20.0], 'Kmag': [12.0]}).to_csv(
        path, index=False)
    window.load_targets(str(path))
    choice = window.target_cols.choice()
    assert choice['shortname'] == 'designation' and choice['dec'] == 'dej2000' and choice['ra'] is None
    assert 'Choose the column for RA' in window.table_status.text() and not window.run_button.isEnabled()
    window.target_cols.boxes['ra'].setCurrentText('alpha_deg')
    assert window.table_status.text() == '1 targets ready.' and window.run_button.isEnabled()
    window.target_list.setText(str(path))
    assert window.kwargs()['target_columns']['ra'] == 'alpha_deg'


def test_charts_appear_as_made(window, tmp_path):
    from matplotlib.figure import Figure
    png = tmp_path / 'chart.png'
    fig = Figure()
    fig.add_subplot().plot([0, 1], [0, 1])
    fig.savefig(png)
    window.add_figure('finder', 'A', fig)
    window.add_figure('finder', 'B', png)  # an Aladin PNG
    assert len(window.items['finder']) == 2 and window.tabs.tabText(2) == 'Finder charts (2)'
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    window.charts.setCurrentRow(1)  # the PNG, drawn in the one canvas
    assert window.chart_view.count() == 1 and window.chart_view.widget(0).findChild(FigureCanvasQTAgg)
    window.charts.setCurrentRow(0)
    canvas = window.chart_view.widget(0).findChild(FigureCanvasQTAgg)
    assert canvas.figure is fig
    canvas.figure.axes[0].set_xlim(0.2, 0.4)  # what the toolbar's zoom does
    canvas.draw()
    assert canvas.figure.axes[0].get_xlim() == (0.2, 0.4)


def test_configure_qt(monkeypatch):
    """A conda qt6.conf must not send PySide6 to another Qt's plugins; explicit user paths win; Wayland
    text-input warnings are silenced without dropping the user's own logging rules."""
    import PySide6

    from targetlist_generator.gui import configure_qt
    monkeypatch.delenv('QT_PLUGIN_PATH', raising=False)
    monkeypatch.setenv('QT_QPA_PLATFORM_PLUGIN_PATH', '/my/plugins')
    monkeypatch.setenv('QT_LOGGING_RULES', 'qt.foo=true')
    configure_qt()
    assert os.environ['QT_PLUGIN_PATH'] == str(Path(PySide6.__file__).parent / 'Qt' / 'plugins')
    assert os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] == '/my/plugins'
    assert os.environ['QT_LOGGING_RULES'] == 'qt.foo=true;qt.qpa.wayland.textinput=false'


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
            assert window.target_list.text().endswith('.vot') and window.target_cols.choice()['mag'] == 'tmassk'
            assert window.from_samp
        finally:
            topcat.disconnect()
    finally:
        window.listen.setChecked(False)
        os.environ.pop('SAMP_HUB', None)
        hub.stop()


def test_matplotlib_uses_pyside6(tmp_path):
    """matplotlib must embed with PySide6 even if another binding is preferred (conda base ships PyQt6)."""
    import subprocess
    import sys
    code = ('from targetlist_generator.gui import figure_page\n'
            'from PySide6.QtWidgets import QApplication\n'
            'from matplotlib.figure import Figure\n'
            'app = QApplication([])\n'
            'figure_page(Figure())\n'
            'from matplotlib.backends import qt_compat\n'
            'print(qt_compat.QT_API)\n')
    env = {**os.environ, 'QT_API': 'pyqt6', 'QT_QPA_PLATFORM': 'offscreen'}
    out = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().endswith('PySide6')


def test_priority_column_picked_up(window, target_csv, tmp_path):
    from targetlist_generator.gui import ColumnMap
    path = tmp_path / 'pri.csv'
    pd.read_csv(target_csv).assign(Prio=2).to_csv(path, index=False)
    window.load_targets(str(path))
    assert window.target_cols.columns()['priority'] == 'prio'  # guessed, and chosen among the other columns
    window.target_cols.boxes['priority'].setCurrentText(ColumnMap.NONE)
    assert window.kwargs()['target_columns']['priority'] is None
    assert 'priority' not in window.std_cols.boxes  # standards have no priority


def test_own_standards_in_window(window, tmp_path):
    assert window.std_cols.choice()['mag'] == 'k'  # bundled list
    path = tmp_path / 'mine.csv'
    pd.DataFrame({'star': ['S1'], 'ra': [10.0], 'dec': [5.0], 'jmag': [7.5], 'kmag': [7.0]}).to_csv(path, index=False)
    window.load_standards_file(str(path))
    assert window.std_cols.choice()['mag'] == 'kmag' and window.std_cols.choice()['shortname'] == 'star'
    window.std_cols.boxes['mag'].setCurrentText('jmag')
    window.target_list.setText('x.csv')
    kw = window.kwargs()
    assert kw['standards_file'] == str(path) and kw['std_band'] == 'jmag' and kw['std_columns']['shortname'] == 'star'
    window.std_file.clear()
    assert window.std_cols.choice()['mag'] == 'k'
