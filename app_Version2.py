import sys
import os
from PyQt5 import QtWidgets, QtCore
from fusion import parse_gtf, build_fusion
import pandas as pd

class FusionApp(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Fusion Transcript Builder")
        self.resize(900, 500)
        layout = QtWidgets.QVBoxLayout(self)

        # file selectors
        h1 = QtWidgets.QHBoxLayout()
        self.fa_edit = QtWidgets.QLineEdit(); btn_fa = QtWidgets.QPushButton("Select FASTA")
        btn_fa.clicked.connect(self.select_fa)
        self.gtf_edit = QtWidgets.QLineEdit(); btn_gtf = QtWidgets.QPushButton("Select GTF")
        btn_gtf.clicked.connect(self.select_gtf)
        h1.addWidget(QtWidgets.QLabel("FASTA:")); h1.addWidget(self.fa_edit); h1.addWidget(btn_fa)
        h1.addSpacing(10)
        h1.addWidget(QtWidgets.QLabel("GTF:")); h1.addWidget(self.gtf_edit); h1.addWidget(btn_gtf)
        layout.addLayout(h1)

        # single input
        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel("Left gene:"),0,0)
        self.left_gene = QtWidgets.QLineEdit(); grid.addWidget(self.left_gene,0,1)
        grid.addWidget(QtWidgets.QLabel("Left breakpoint (chr:pos):"),0,2)
        self.left_bp = QtWidgets.QLineEdit(); grid.addWidget(self.left_bp,0,3)

        grid.addWidget(QtWidgets.QLabel("Right gene:"),1,0)
        self.right_gene = QtWidgets.QLineEdit(); grid.addWidget(self.right_gene,1,1)
        grid.addWidget(QtWidgets.QLabel("Right breakpoint (chr:pos):"),1,2)
        self.right_bp = QtWidgets.QLineEdit(); grid.addWidget(self.right_bp,1,3)
        layout.addLayout(grid)

        # batch input
        h2 = QtWidgets.QHBoxLayout()
        self.csv_edit = QtWidgets.QLineEdit(); btn_csv = QtWidgets.QPushButton("Batch CSV")
        btn_csv.clicked.connect(self.select_csv)
        h2.addWidget(QtWidgets.QLabel("Batch CSV:")); h2.addWidget(self.csv_edit); h2.addWidget(btn_csv)
        layout.addLayout(h2)

        # run & output
        h3 = QtWidgets.QHBoxLayout()
        btn_run = QtWidgets.QPushButton("Run")
        btn_run.clicked.connect(self.run)
        self.save_btn = QtWidgets.QPushButton("Export TSV")
        self.save_btn.clicked.connect(self.export_tsv)
        self.save_btn.setEnabled(False)
        h3.addWidget(btn_run); h3.addWidget(self.save_btn)
        layout.addLayout(h3)

        # results table
        self.table = QtWidgets.QTableWidget()
        layout.addWidget(self.table)

        self.gtf_dict = None
        self.fasta_path = None
        self.results = []

    def select_fa(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select FASTA", filter="FASTA (*.fa *.fasta *.fa.gz *.fna)")
        if p:
            self.fa_edit.setText(p); self.fasta_path = p

    def select_gtf(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select GTF", filter="GTF (*.gtf *.gtf.gz)")
        if p:
            self.gtf_edit.setText(p)
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            try:
                self.gtf_dict = parse_gtf(p)
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, "GTF parse error", str(e))
                self.gtf_dict = None
            finally:
                QtWidgets.QApplication.restoreOverrideCursor()

    def select_csv(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select CSV", filter="CSV (*.csv)")
        if p:
            self.csv_edit.setText(p)

    def run(self):
        if not self.fasta_path or not self.gtf_dict:
            QtWidgets.QMessageBox.warning(self, "Missing files", "Please select FASTA and GTF first.")
            return
        rows = []
        if self.csv_edit.text():
            try:
                df = pd.read_csv(self.csv_edit.text())
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, "CSV error", str(e))
                return
            # 尝试容错列名大小写
            cols_lower = {c.lower(): c for c in df.columns}
            required = ['leftgene','rightgene','leftbp','rightbp']
            if not all(k in cols_lower for k in required):
                QtWidgets.QMessageBox.critical(self, "CSV format", "CSV must contain columns: leftgene,rightgene,leftbp,rightbp")
                return
            for _, r in df.iterrows():
                rows.append((r[cols_lower['leftgene']], r[cols_lower['rightgene']], r[cols_lower['leftbp']], r[cols_lower['rightbp']]))
        else:
            rows.append((self.left_gene.text().strip(), self.right_gene.text().strip(), self.left_bp.text().strip(), self.right_bp.text().strip()))
        self.results = []
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            for lg, rg, lb, rb in rows:
                try:
                    res = build_fusion(self.fasta_path, self.gtf_dict, lg, rg, lb, rb)
                    self.results.append(res)
                except Exception as e:
                    self.results.append({'left_gene':lg,'right_gene':rg,'left_tx':'','right_tx':'','fusion_nt':'','fusion_protein':'','notes':[str(e)]})
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.populate_table()
        self.save_btn.setEnabled(True)

    def populate_table(self):
        cols = ['left_gene','right_gene','left_tx','right_tx','fusion_nt','fusion_protein','notes']
        self.table.setColumnCount(len(cols)); self.table.setRowCount(len(self.results))
        self.table.setHorizontalHeaderLabels(cols)
        for i,res in enumerate(self.results):
            for j,c in enumerate(cols):
                v = res.get(c, '')
                if isinstance(v, list): v = '; '.join(v)
                it = QtWidgets.QTableWidgetItem(str(v))
                self.table.setItem(i,j,it)
        self.table.resizeColumnsToContents()

    def export_tsv(self):
        p, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save TSV", filter="TSV (*.tsv)")
        if not p: return
        import csv
        cols = ['left_gene','right_gene','left_tx','right_tx','fusion_nt','fusion_protein','notes']
        with open(p, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
            writer.writeheader()
            for r in self.results:
                rr = {k: ('; '.join(r[k]) if isinstance(r.get(k), list) else r.get(k,'')) for k in cols}
                writer.writerow(rr)
        QtWidgets.QMessageBox.information(self, "Saved", f"Saved to {p}")

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = FusionApp()
    w.show()
    sys.exit(app.exec_())