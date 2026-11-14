"""Graphical interface: load a table (file or SAMP), set times, night and observatory, plan, save."""
import os
import sys
import tempfile
from pathlib import Path
from urllib.request import urlretrieve

import PySide6
from PySide6.QtCore import QDate, QObject, Qt, QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QCompleter,
    QDateEdit,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from . import finder, formats, io, sky, standards
from .pipeline import PLAN_COLUMNS, make_targetlist

CUSTOM = 'Custom (latitude/longitude below)'
DEFAULT_SITE = 'NASA Infrared Telescope Facility'
TABLES = 'Tables (*.csv *.fits *.fit *.ecsv *.vot *.xml *.txt);;All files (*)'
SAMP_TYPES = {'table.load.votable': '.vot', 'table.load.fits': '.fits'}


def spin(lo: float, hi: float, value: float, decimals: int = 2, step: float = 1.0) -> QDoubleSpinBox:
    box = QDoubleSpinBox(decimals=decimals, minimum=lo, maximum=hi, singleStep=step)
    box.setValue(value)
    return box


def form(group: QGroupBox, rows) -> QGroupBox:
    layout = QFormLayout(group)
    for label, widget in rows:
        layout.addRow(label, widget)
    return group


def pair(*widgets) -> QHBoxLayout:
    row = QHBoxLayout()
    for widget in widgets:
        row.addWidget(widget)
    return row


class Worker(QObject):
    """Runs the pipeline off the GUI thread."""
    log = Signal(str)
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, kwargs: dict):
        super().__init__()
        self.kwargs = kwargs

    def run(self):
        try:
            self.done.emit(make_targetlist(log=self.log.emit, **self.kwargs))
        except Exception as err:  # thread boundary: every failure is shown to the user, never swallowed
            self.failed.emit(f'{type(err).__name__}: {err}')


class Samp(QObject):
    """Receives tables sent from TOPCAT (or any SAMP application); starts its own hub when none is running."""
    table = Signal(str, str)  # local path, table name
    failed = Signal(str)

    def __init__(self):
        super().__init__()
        self.client = self.hub = None
        self.folder = tempfile.TemporaryDirectory(prefix='targetlist_samp_')

    def start(self):
        from pyvo.samp import SAMPHubError, SAMPHubServer, SAMPIntegratedClient
        self.client = SAMPIntegratedClient(name='Target List Generator',
                                           description='Plans a night from a target table')
        try:
            self.client.connect()
        except SAMPHubError:  # no hub yet: be the hub, so TOPCAT started later joins us
            self.hub = SAMPHubServer(web_profile=False)
            self.hub.start()
            self.client.connect()
        for mtype in SAMP_TYPES:
            self.client.bind_receive_notification(mtype, self.notified)
            self.client.bind_receive_call(mtype, self.called)

    def stop(self):
        if self.client is not None and self.client.is_connected:
            self.client.disconnect()
        if self.hub is not None:
            self.hub.stop()
        self.client = self.hub = None

    def fetch(self, mtype: str, params: dict):
        name = params.get('name') or 'samp_table'
        path = Path(self.folder.name) / f'{len(list(Path(self.folder.name).iterdir()))}{SAMP_TYPES[mtype]}'
        try:
            urlretrieve(params['url'], path)
        except (OSError, KeyError) as err:
            self.failed.emit(f'Could not fetch table {name!r} sent by SAMP: {err}')
            return
        self.table.emit(str(path), name)

    def notified(self, private_key, sender_id, mtype, params, extra):
        self.fetch(mtype, params)

    def called(self, private_key, sender_id, msg_id, mtype, params, extra):
        self.fetch(mtype, params)
        self.client.reply(msg_id, {'samp.status': 'samp.ok', 'samp.result': {}})


class MainWindow(QMainWindow):
    finished = Signal()  # after each run, success or not

    def __init__(self):
        super().__init__()
        self.setWindowTitle('Target List Generator')
        self.thread = self.worker = None
        self.samp = Samp()
        self.samp.table.connect(self.samp_table)
        self.samp.failed.connect(lambda msg: QMessageBox.warning(self, 'SAMP', msg))

        # 1. table
        self.target_list = QLineEdit(placeholderText='CSV, FITS, ECSV, VOTable or IRTF .txt, or send from TOPCAT')
        self.listen = QCheckBox('Listen for TOPCAT (SAMP)')
        self.listen.toggled.connect(self.toggle_samp)
        self.mag_col = QComboBox(editable=True)
        table = QGroupBox('1. Target table')
        table_form = QFormLayout(table)
        table_form.addRow('File', self.with_browse(self.target_list, self.browse_targets))
        table_form.addRow('', self.listen)
        table_form.addRow('Magnitude column', self.mag_col)

        # 2. times
        self.time_same, self.time_column = QRadioButton('Same for all [min]'), QRadioButton('From column')
        self.time_same.setChecked(True)
        self.target_time = spin(1, 600, 30, 0, 5)
        self.time_col = QComboBox()
        self.std_time = spin(0, 120, 10, 0, 5)
        self.cal_time = spin(0, 120, 5, 0, 5)
        times = form(QGroupBox('2. Time per target block'), (
            (self.time_same, self.target_time), (self.time_column, self.time_col),
            ('Standard [min]', self.std_time), ('Calibrations [min]', self.cal_time)))

        # 3. night
        self.night = QGroupBox('3. Night and observatory', checkable=True)
        self.date = QDateEdit(QDate.currentDate(), calendarPopup=True, displayFormat='yyyy-MM-dd')
        self.site = QComboBox(editable=True, insertPolicy=QComboBox.NoInsert)
        sites = sky.observatories()
        self.site.addItem(CUSTOM, None)
        for i, row in sites.iterrows():
            aliases = row['aliases'].replace(';', ', ')
            self.site.addItem(row['name'] + (f'  [{aliases}]' if aliases else ''), i)
        completer = QCompleter(self.site.model(), self.site)
        completer.setFilterMode(Qt.MatchContains)
        completer.setCaseSensitivity(Qt.CaseInsensitive)
        self.site.setCompleter(completer)
        self.lat, self.lon = spin(-90, 90, 0, 4), spin(-180, 180, 0, 4)
        self.site.currentIndexChanged.connect(self.site_changed)
        for box in (self.lat, self.lon):
            box.valueChanged.connect(self.coords_typed)
        self.twilight = spin(-18, 0, -12, 1)
        self.max_airmass = spin(1, 5, 2, 2, 0.1)
        self.min_time = spin(0, 720, 30, 0, 10)
        self.moon = QComboBox()
        self.moon.addItems(('auto (10-40 deg with phase)', 'off', 'custom'))
        self.moon_sep = spin(0, 180, 30, 0, 5)
        self.moon_sep.setEnabled(False)
        self.moon.currentTextChanged.connect(lambda t: self.moon_sep.setEnabled(t == 'custom'))
        form(self.night, (('Evening of', self.date), ('Observatory (type to search)', self.site),
                          ('Latitude [deg]', self.lat), ('Longitude [deg, east +]', self.lon),
                          ('Twilight Sun altitude [deg]', self.twilight), ('Max airmass', self.max_airmass),
                          ('Min time below it [min]', self.min_time),
                          ('Moon separation', pair(self.moon, self.moon_sep))))
        self.site.setCurrentIndex(self.site.findText(DEFAULT_SITE, Qt.MatchStartsWith))

        # 4. standards
        self.stds = QGroupBox('4. Telluric standards', checkable=True)
        self.std_band, self.std_spt = QComboBox(), QComboBox()
        self.std_band.addItems(standards.BANDS)
        self.std_band.setCurrentText('K')
        self.std_spt.addItems(('any', 'A0V', 'A1V', 'B9V'))
        self.std_tol = spin(0.01, 1, 0.1, 2, 0.05)
        self.std_sep = spin(1, 90, 20, 0, 5)
        form(self.stds, (('Band', self.std_band), ('Spectral type', self.std_spt),
                         ('Airmass tolerance', self.std_tol), ('Max distance [deg]', self.std_sep)))

        # 5. optional RA/Dec box
        self.box = QGroupBox('5. RA/Dec box (optional)', checkable=True, checked=False)
        self.ra_lo, self.ra_hi = QLineEdit('20h'), QLineEdit('4h')
        self.dec_lo, self.dec_hi = QLineEdit('-30'), QLineEdit('60')
        form(self.box, (('RA from / to (deg or 20h)', pair(self.ra_lo, self.ra_hi)),
                        ('Dec from / to (deg)', pair(self.dec_lo, self.dec_hi))))

        # 6. outputs
        outputs = QGroupBox('6. Outputs')
        grid = QGridLayout()
        self.format_boxes = {}
        for n, (key, (label, _)) in enumerate(formats.FORMATS.items()):
            self.format_boxes[key] = QCheckBox(label, checked=key in formats.DEFAULT)
            grid.addWidget(self.format_boxes[key], n // 2, n % 2)
        self.finders = QCheckBox('Finder charts')
        self.survey = QComboBox(editable=True)
        self.survey.addItems(finder.SURVEYS)
        self.fov = spin(0.2, 60, 2, 1, 0.5)
        self.mirror = QLineEdit(placeholderText='optional hips2fits mirror URL')
        self.out_dir = QLineEdit(placeholderText='asked when you press Run')
        out_form = QFormLayout(outputs)
        out_form.addRow(grid)
        out_form.addRow(self.finders, pair(self.survey, QLabel('width [arcmin]'), self.fov))
        out_form.addRow('Chart mirror', self.mirror)
        out_form.addRow('Save in', self.with_browse(self.out_dir, self.browse_out))

        self.run_button = QPushButton('Plan the night')
        self.run_button.clicked.connect(self.run)
        open_button = QPushButton('Open output folder')
        open_button.clicked.connect(self.open_output)

        inputs = QWidget()
        column = QVBoxLayout(inputs)
        for widget in (table, times, self.night, self.stds, self.box, outputs):
            column.addWidget(widget)
        column.addLayout(pair(self.run_button, open_button))
        column.addStretch()
        scroll = QScrollArea(widgetResizable=True)
        scroll.setWidget(inputs)

        self.log = QPlainTextEdit(readOnly=True)
        self.table = QTableWidget()
        self.plots = QComboBox()
        self.plot = QLabel('Give a night to see visibility plots', alignment=Qt.AlignCenter)
        self.plots.currentIndexChanged.connect(self.show_plot)
        self.plot_scroll = QScrollArea(widgetResizable=True)
        self.plot_scroll.setWidget(self.plot)
        plot_view = QWidget()
        plot_layout = QVBoxLayout(plot_view)
        plot_layout.addWidget(self.plots)
        plot_layout.addWidget(self.plot_scroll)
        self.charts = QListWidget()
        self.chart = QLabel(alignment=Qt.AlignCenter)
        self.charts.currentTextChanged.connect(self.show_chart)
        chart_view = QSplitter()
        chart_view.addWidget(self.charts)
        chart_view.addWidget(self.chart)
        self.tabs = QTabWidget()
        self.tabs.addTab(self.table, 'Plan')
        self.tabs.addTab(plot_view, 'Visibility')
        self.tabs.addTab(chart_view, 'Finder charts')
        self.tabs.addTab(self.log, 'Log')
        self.tabs.currentChanged.connect(lambda _: self.show_plot(self.plots.currentIndex()))  # rescale to fit

        main = QSplitter()
        main.addWidget(scroll)
        main.addWidget(self.tabs)
        main.setSizes([460, 760])
        self.setCentralWidget(main)
        self.resize(1250, 820)
        self.plot_paths: list[Path] = []
        self.chart_paths: dict[str, Path] = {}

    def with_browse(self, line: QLineEdit, slot) -> QHBoxLayout:
        button = QPushButton('Browse…')
        button.clicked.connect(slot)
        return pair(line, button)

    def browse_targets(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Target table', str(Path.home()), TABLES)
        if path:
            self.load_targets(path)

    def browse_out(self):
        path = QFileDialog.getExistingDirectory(self, 'Where should the lists and charts be saved?',
                                                self.out_dir.text() or str(Path.home()))
        if path:
            self.out_dir.setText(path)

    def load_targets(self, path: str, label: str | None = None):
        """Set the table and offer its numeric columns as magnitude and time columns."""
        try:
            df = io.read_table(path)
        except (OSError, ValueError) as err:
            QMessageBox.critical(self, 'Cannot read table', str(err))
            return
        self.target_list.setText(path)
        numeric = list(df.select_dtypes('number').columns.drop(['ra', 'dec'], errors='ignore'))
        for box in (self.mag_col, self.time_col):
            box.clear()
            box.addItems(numeric)
        self.mag_col.setCurrentText('tmassk' if 'tmassk' in df.columns else self.mag_col.itemText(0))
        self.log.appendPlainText(f'{label or path}: {len(df)} rows, columns {list(df.columns)}')

    def toggle_samp(self, on: bool):
        try:
            self.samp.start() if on else self.samp.stop()
        except Exception as err:  # hub/network failures are shown, and the toggle reset
            self.listen.blockSignals(True)
            self.listen.setChecked(False)
            self.listen.blockSignals(False)
            QMessageBox.critical(self, 'SAMP', f'Could not connect: {err}')
            return
        self.log.appendPlainText('Listening for tables from TOPCAT (SAMP)' if on else 'Stopped listening (SAMP)')

    def samp_table(self, path: str, name: str):
        self.load_targets(path, label=f'From SAMP: {name}')
        self.raise_()

    def site_changed(self, index: int):
        row = self.site.itemData(index)
        if row is None:
            return
        site = sky.observatories().loc[row]
        for box, value in ((self.lat, site['lat']), (self.lon, site['lon'])):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)

    def coords_typed(self):
        self.site.blockSignals(True)
        self.site.setCurrentIndex(0)  # typed coordinates: custom site
        self.site.blockSignals(False)

    def kwargs(self) -> dict:
        if not self.target_list.text():
            raise ValueError('Choose a target table first (or send one from TOPCAT)')
        keys = [k for k, box in self.format_boxes.items() if box.isChecked()]
        if not keys:
            raise ValueError('Tick at least one output format')
        mirror = self.mirror.text().strip()
        kw = {'target_list': self.target_list.text(), 'out_dir': self.out_dir.text(),
              'mag_col': self.mag_col.currentText(), 'target_time': self.target_time.value(),
              'target_time_col': self.time_col.currentText() if self.time_column.isChecked() else None,
              'cal_time': self.cal_time.value(), 'std_time': self.std_time.value(),
              'with_standards': self.stds.isChecked(), 'std_band': self.std_band.currentText(),
              'std_spt': None if self.std_spt.currentText() == 'any' else self.std_spt.currentText(),
              'std_airmass_tol': self.std_tol.value(), 'std_max_sep': self.std_sep.value(), 'outputs': keys,
              'finders': self.finders.isChecked(), 'survey': self.survey.currentText(), 'fov': self.fov.value(),
              'urls': ((mirror,) if mirror else ()) + finder.HIPS2FITS}
        if self.time_column.isChecked() and not self.time_col.currentText():
            raise ValueError('Choose the column holding each target\'s time')
        if self.night.isChecked():
            moon = self.moon.currentText()
            kw |= {'date': self.date.date().toString('yyyy-MM-dd'), 'twilight': self.twilight.value(),
                   'max_airmass': self.max_airmass.value(), 'min_time': self.min_time.value(),
                   'min_moon_sep': self.moon_sep.value() if moon == 'custom' else moon.split()[0]}
            row = self.site.currentData()
            if row is None:
                kw |= {'lat': self.lat.value(), 'lon': self.lon.value()}
            else:
                kw['site'] = sky.observatories().loc[row, 'name']
        if self.box.isChecked():
            kw |= {'ra_range': (self.ra_lo.text(), self.ra_hi.text()),
                   'dec_range': (self.dec_lo.text(), self.dec_hi.text())}
        return kw

    def run(self):
        try:
            kwargs = self.kwargs()
        except ValueError as err:
            QMessageBox.warning(self, 'Check inputs', str(err))
            return
        if not kwargs['out_dir']:
            self.browse_out()
            if not self.out_dir.text():
                return
            kwargs['out_dir'] = self.out_dir.text()
        self.run_button.setEnabled(False)
        self.tabs.setCurrentWidget(self.log)
        self.thread, self.worker = QThread(self), Worker(kwargs)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.log.connect(self.log.appendPlainText)
        self.worker.done.connect(self.on_done)
        self.worker.failed.connect(self.on_failed)
        self.thread.start()

    def end_run(self):
        self.thread.quit()
        self.thread.wait()
        self.run_button.setEnabled(True)
        self.finished.emit()

    def on_done(self, result):
        rows, paths = result
        if 'status' in rows.columns:  # planned night
            shown = rows[[c for c in PLAN_COLUMNS if c in rows.columns and c not in ('ra', 'dec')]]
        else:
            shown = io.formatted(rows, self.mag_col.currentText().lower())
        self.table.setColumnCount(len(shown.columns))
        self.table.setRowCount(len(shown))
        self.table.setHorizontalHeaderLabels(list(shown.columns))
        for i, row in enumerate(shown.itertuples(index=False)):
            for j, value in enumerate(row):
                self.table.setItem(i, j, QTableWidgetItem('' if value is None else str(value)))
        self.table.resizeColumnsToContents()
        self.plot_paths = [p for p in paths if '_visibility_' in p.name]
        self.plots.clear()
        self.plots.addItems([p.stem for p in self.plot_paths])
        self.chart_paths = {p.stem: p for p in paths if p.parent.name == 'fcharts'}
        self.charts.clear()
        self.charts.addItems(list(self.chart_paths))
        self.tabs.setCurrentWidget(self.table)
        self.end_run()

    def on_failed(self, message: str):
        self.log.appendPlainText(message)
        self.end_run()
        QMessageBox.critical(self, 'Run failed', message)

    def show_plot(self, index: int):
        if 0 <= index < len(self.plot_paths):
            width = self.plot_scroll.viewport().width() - 4  # fit the tab; scroll for the rest
            self.plot.setPixmap(QPixmap(str(self.plot_paths[index])).scaledToWidth(width, Qt.SmoothTransformation))

    def show_chart(self, name: str):
        if name in self.chart_paths:
            self.chart.setPixmap(QPixmap(str(self.chart_paths[name])).scaled(
                self.chart.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation))

    def open_output(self):
        if self.out_dir.text():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(self.out_dir.text()).resolve())))

    def closeEvent(self, event):
        self.samp.stop()
        super().closeEvent(event)


def use_bundled_qt_plugins() -> None:
    """Point Qt at the plugins shipped with PySide6 (pip wheels, frozen apps). Otherwise a qt6.conf next to the
    Python executable, e.g. from conda's qtbase package, sends Qt to another Qt version's platform plugins and the
    window aborts with 'no Qt platform plugin could be initialized'. Paths the user set explicitly win."""
    plugins = Path(PySide6.__file__).parent / 'Qt' / 'plugins'
    if (plugins / 'platforms').is_dir():  # system-packaged PySide6 uses the system Qt's own plugins instead
        os.environ.setdefault('QT_PLUGIN_PATH', str(plugins))
        os.environ.setdefault('QT_QPA_PLATFORM_PLUGIN_PATH', str(plugins / 'platforms'))


def main():
    use_bundled_qt_plugins()
    app = QApplication.instance() or QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
