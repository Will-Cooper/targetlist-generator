"""Graphical interface: pick a target list, a night and a site, press Run."""
import sys
from pathlib import Path

from PySide6.QtCore import QDate, QObject, Qt, QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDateEdit,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from . import finder, io, standards
from .pipeline import make_targetlist

CUSTOM = 'Custom latitude/longitude'
SITES = (CUSTOM, 'irtf', 'lapalma', 'paranal', 'lasilla', 'Cerro Pachon', 'gemini_north', 'keck', 'Kitt Peak')
TABLES = 'Tables (*.csv *.fits *.fit *.ecsv *.vot *.xml *.txt);;All files (*)'


def spin(lo: float, hi: float, value: float, decimals: int = 2, step: float = 1.0) -> QDoubleSpinBox:
    box = QDoubleSpinBox(decimals=decimals, minimum=lo, maximum=hi, singleStep=step)
    box.setValue(value)
    return box


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


class MainWindow(QMainWindow):
    finished = Signal()  # after each run, success or not

    def __init__(self):
        super().__init__()
        self.setWindowTitle('Target List Generator')
        self.thread = self.worker = None

        self.target_list = QLineEdit(placeholderText='CSV, FITS, ECSV, VOTable or IRTF .txt')
        self.mag_col = QComboBox(editable=True)
        self.out_dir = QLineEdit(str(Path.cwd() / 'targetlist_out'))

        self.night = QGroupBox('Night visibility', checkable=True)
        self.date = QDateEdit(QDate.currentDate(), calendarPopup=True, displayFormat='yyyy-MM-dd')
        self.site = QComboBox(editable=True)
        self.site.addItems(SITES)
        self.lat, self.lon = spin(-90, 90, 0, 4), spin(-180, 180, 0, 4)
        self.twilight = spin(-18, 0, -12, 1)
        self.max_airmass = spin(1, 5, 2, 2, 0.1)
        self.min_time = spin(0, 720, 30, 0, 10)
        self.min_moon_sep = spin(0, 180, 0, 0, 5)
        self.order = QComboBox()
        self.order.addItems(('time', 'input'))
        self.site.currentTextChanged.connect(self.site_changed)
        self.site.setCurrentText('irtf')
        night = QFormLayout(self.night)
        for label, widget in (('Evening of', self.date), ('Site', self.site), ('Latitude [deg]', self.lat),
                              ('Longitude [deg, east +]', self.lon), ('Twilight Sun alt. [deg]', self.twilight),
                              ('Max airmass', self.max_airmass), ('Min time below it [min]', self.min_time),
                              ('Min Moon separation [deg]', self.min_moon_sep),
                              ('List order (best time / input)', self.order)):
            night.addRow(label, widget)

        self.box = QGroupBox('RA/Dec box', checkable=True, checked=False)
        self.ra_lo, self.ra_hi = QLineEdit('20h'), QLineEdit('4h')
        self.dec_lo, self.dec_hi = QLineEdit('-30'), QLineEdit('60')
        box = QFormLayout(self.box)
        box.addRow('RA from / to (deg or 20h)', self.pair(self.ra_lo, self.ra_hi))
        box.addRow('Dec from / to (deg)', self.pair(self.dec_lo, self.dec_hi))

        self.stds = QGroupBox('Telluric standards', checkable=True)
        self.std_band, self.std_spt = QComboBox(), QComboBox()
        self.std_band.addItems(standards.BANDS)
        self.std_band.setCurrentText('K')
        self.std_spt.addItems(('any', 'A0V', 'A1V', 'B9V'))
        self.std_offset = spin(-120, 120, 30, 0, 5)
        self.std_airmass_tol = spin(0.01, 1, 0.1, 2, 0.05)
        stds = QFormLayout(self.stds)
        for label, widget in (('Band', self.std_band), ('Spectral type', self.std_spt),
                              ('Observed after target [min]', self.std_offset),
                              ('Airmass tolerance (with a night)', self.std_airmass_tol)):
            stds.addRow(label, widget)

        self.finders = QGroupBox('Finder charts', checkable=True, checked=False)
        self.survey = QComboBox(editable=True)
        self.survey.addItems(finder.SURVEYS)
        self.fov = spin(0.2, 60, 2, 1, 0.5)
        self.mirror = QLineEdit(placeholderText='optional hips2fits mirror URL')
        fcs = QFormLayout(self.finders)
        for label, widget in (('Survey', self.survey), ('Width [arcmin]', self.fov), ('Mirror', self.mirror)):
            fcs.addRow(label, widget)

        self.run_button = QPushButton('Run')
        self.run_button.clicked.connect(self.run)
        open_button = QPushButton('Open output folder')
        open_button.clicked.connect(self.open_output)

        inputs = QWidget()
        form = QVBoxLayout(inputs)
        files = QFormLayout()
        files.addRow('Target list', self.with_browse(self.target_list, self.browse_targets))
        files.addRow('Magnitude column', self.mag_col)
        files.addRow('Output folder', self.with_browse(self.out_dir, self.browse_out))
        form.addLayout(files)
        for widget in (self.night, self.box, self.stds, self.finders):
            form.addWidget(widget)
        form.addLayout(self.pair(self.run_button, open_button))
        form.addStretch()
        scroll = QScrollArea(widgetResizable=True)
        scroll.setWidget(inputs)

        self.log = QPlainTextEdit(readOnly=True)
        self.table = QTableWidget()
        self.charts = QListWidget()
        self.chart = QLabel(alignment=Qt.AlignCenter)
        self.charts.currentTextChanged.connect(self.show_chart)
        chart_view = QSplitter()
        chart_view.addWidget(self.charts)
        chart_view.addWidget(self.chart)
        self.tabs = QTabWidget()
        self.tabs.addTab(self.log, 'Log')
        self.tabs.addTab(self.table, 'Targets')
        self.tabs.addTab(chart_view, 'Finder charts')
        self.airmass = QLabel('Give a night to see the airmass plot', alignment=Qt.AlignCenter)
        self.tabs.addTab(self.airmass, 'Airmass plot')

        main = QSplitter()
        main.addWidget(scroll)
        main.addWidget(self.tabs)
        main.setSizes([420, 700])
        self.setCentralWidget(main)
        self.resize(1150, 750)
        self.chart_paths: dict[str, Path] = {}

    @staticmethod
    def pair(*widgets) -> QHBoxLayout:
        row = QHBoxLayout()
        for widget in widgets:
            row.addWidget(widget)
        return row

    def with_browse(self, line: QLineEdit, slot) -> QHBoxLayout:
        button = QPushButton('Browse…')
        button.clicked.connect(slot)
        return self.pair(line, button)

    def browse_targets(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Target list', str(Path.cwd()), TABLES)
        if path:
            self.load_targets(path)

    def browse_out(self):
        path = QFileDialog.getExistingDirectory(self, 'Output folder', self.out_dir.text())
        if path:
            self.out_dir.setText(path)

    def load_targets(self, path: str):
        """Set the list and offer its numeric columns as the magnitude column."""
        try:
            df = io.read_table(path)
        except (OSError, ValueError) as err:
            QMessageBox.critical(self, 'Cannot read table', str(err))
            return
        self.target_list.setText(path)
        self.mag_col.clear()
        self.mag_col.addItems(list(df.select_dtypes('number').columns.drop(['ra', 'dec'], errors='ignore')))
        self.mag_col.setCurrentText('tmassk' if 'tmassk' in df.columns else self.mag_col.itemText(0))
        self.log.appendPlainText(f'{path}: {len(df)} rows, columns {list(df.columns)}')

    def site_changed(self, text: str):
        for widget in (self.lat, self.lon):
            widget.setEnabled(text == CUSTOM)

    def kwargs(self) -> dict:
        if not self.target_list.text():
            raise ValueError('Choose a target list first')
        kw = {'target_list': self.target_list.text(), 'out_dir': self.out_dir.text(),
              'mag_col': self.mag_col.currentText(), 'with_standards': self.stds.isChecked(),
              'std_band': self.std_band.currentText(), 'std_offset': self.std_offset.value(),
              'std_airmass_tol': self.std_airmass_tol.value(), 'order': self.order.currentText(),
              'std_spt': None if self.std_spt.currentText() == 'any' else self.std_spt.currentText(),
              'finders': self.finders.isChecked(), 'survey': self.survey.currentText(), 'fov': self.fov.value(),
              'urls': ((self.mirror.text().strip(),) if self.mirror.text().strip() else ()) + finder.HIPS2FITS}
        if self.night.isChecked():
            kw |= {'date': self.date.date().toString('yyyy-MM-dd'), 'twilight': self.twilight.value(),
                   'max_airmass': self.max_airmass.value(), 'min_time': self.min_time.value(),
                   'min_moon_sep': self.min_moon_sep.value()}
            if self.site.currentText() == CUSTOM:
                kw |= {'lat': self.lat.value(), 'lon': self.lon.value()}
            else:
                kw['site'] = self.site.currentText()
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
        kept, paths = result
        rows = io.formatted(kept, self.mag_col.currentText().lower())
        for col in ('best_utc', 'min_airmass', 'moon_sep', 'standard'):  # present when a night/standards given
            if col in kept.columns:
                rows[col] = kept[col].to_numpy()
        self.table.setColumnCount(len(rows.columns))
        self.table.setRowCount(len(rows))
        self.table.setHorizontalHeaderLabels(list(rows.columns))
        for i, row in enumerate(rows.itertuples(index=False)):
            for j, value in enumerate(row):
                self.table.setItem(i, j, QTableWidgetItem(str(value)))
        plots = [p for p in paths if p.name.endswith('_airmass.png')]
        if plots:
            self.airmass.setPixmap(QPixmap(str(plots[0])))
        self.chart_paths = {p.stem: p for p in paths if p.suffix == '.png' and p not in plots}
        self.charts.clear()
        self.charts.addItems(list(self.chart_paths))
        self.tabs.setCurrentWidget(self.table)
        self.end_run()

    def on_failed(self, message: str):
        self.log.appendPlainText(message)
        self.end_run()
        QMessageBox.critical(self, 'Run failed', message)

    def show_chart(self, name: str):
        if name in self.chart_paths:
            self.chart.setPixmap(QPixmap(str(self.chart_paths[name])).scaled(
                self.chart.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation))

    def open_output(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(self.out_dir.text()).resolve())))


def main():
    app = QApplication.instance() or QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
