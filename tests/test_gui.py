import pytest

pytest.importorskip('PySide6')
pytestmark = pytest.mark.gui


def test_gui_run(qtbot, target_csv, tmp_path):
    from PySide6.QtCore import QDate

    from targetlist_generator.gui import CUSTOM, MainWindow
    window = MainWindow()
    qtbot.addWidget(window)
    window.load_targets(str(target_csv))
    assert window.mag_col.currentText() == 'tmassk'
    window.site.setCurrentText(CUSTOM)
    window.lat.setValue(19.83)
    window.lon.setValue(-155.47)
    window.date.setDate(QDate(2026, 11, 14))
    window.out_dir.setText(str(tmp_path))
    with qtbot.waitSignal(window.finished, timeout=60_000):
        window.run_button.click()
    assert 0 < window.table.rowCount() < 25
    assert (tmp_path / 'targets_irtf.txt').exists()
    headers = [window.table.horizontalHeaderItem(j).text() for j in range(window.table.columnCount())]
    assert {'best_utc', 'min_airmass', 'moon_sep', 'standard'} <= set(headers)
    assert not window.airmass.pixmap().isNull()


def test_gui_needs_target_list(qtbot, monkeypatch):
    from PySide6.QtWidgets import QMessageBox

    from targetlist_generator.gui import MainWindow
    shown = []
    monkeypatch.setattr(QMessageBox, 'warning', lambda *args: shown.append(args[-1]))
    window = MainWindow()
    qtbot.addWidget(window)
    window.run_button.click()
    assert shown == ['Choose a target list first']
