#-------------------------------------------------------------------------------
# Contains handlers for GenomicsDB search end-points
#
#
#-------------------------------------------------------------------------------
import os
import time
import logging
from flask import render_template, flash, redirect, session, url_for, request, g
from flask import make_response, send_from_directory, send_file
from werkzeug.utils import secure_filename
from app import app
from app import ds
from forms import SearchForm, SearchRangeForm, UploadForm, DownloadForm
from datetime import datetime
from utils import allowed_file

@app.errorhandler(404)
def internal_error(error):
  return render_template('404.html'), 404

#@app.errorhandler(500)
#def internal_error(error):
#  return render_template('500.html'), 500

@app.route('/', methods = ['GET', 'POST'])
@app.route('/index', methods = ['GET', 'POST'])
@app.route('/index/<rsid>', methods = ['GET', 'POST'])
@app.route('/index/<rsid>/<threshold>', methods = ['GET', 'POST'])
@app.route('/index/<rsid>/<threshold>/<procopt>', methods = ['GET', 'POST'])
def index(rsid=None, threshold=0.9, procopt=None):
  """
  Home page - search by single rsid, allow probability threshold adjustment
  """
  logging.info("request method - %s", str(request.method))
  form = SearchForm()
  msg = ""
  if form.validate_on_submit():
    logging.info("validate form procopt - %s", procopt)
    if "dbtn" in request.form and procopt != "download":
      procopt = 'download'
      return redirect(url_for('index', rsid=form.rs.data, threshold=form.threshold.data, procopt=procopt), code=307)
    if procopt != "download":
      procopt = 'display'
      return redirect(url_for('index', rsid=form.rs.data, threshold=form.threshold.data, procopt=procopt))
  variant_data = []
  threshold = float(threshold)
  if rsid != None:
    logging.info("RSID-present %s", rsid)
    if procopt == 'download':
      logging.info("Initiate download")
      logging.info("request - %s", str(request.form))
      download_list = request.form.getlist("include")
      logging.info("download_list - %s", str(download_list))
      (sample_return_data, snp_return_data, msg) = ds.get_rslist_data([rsid], threshold, download_list)
      zipfilename = rsid + "_single_snp_results.zip"
      zipFilePath = ds.make_zipfile(sample_return_data, snp_return_data, app.config['UPLOAD_DIR'], zipfilename)
      return send_file(zipFilePath, as_attachment=True)

    (variant_data, msg) = ds.get_variant_summary_probs(rsid, threshold)
    form.rs.data = rsid
    form.threshold.data = threshold
  logging.info("Drop-through")
  return render_template('index.html',
    title = 'Home',
    form = form,
    variant_data = variant_data,
    msg = msg,
    db_name = ds.get_db_name())

@app.route('/range', methods = ['GET', 'POST'])
@app.route('/range/<chromosome>/<start>/<end>/<procopt>', methods = ['GET', 'POST'])
def range(chromosome = None, start = None, end=None, procopt=None):
  """
  Search by range within chromosome
  """
  form = SearchRangeForm(request.form)
  msg = ""
  threshold = 0.9
  if form.validate_on_submit():
    logging.info("range - form validated")
    #threshold = form.threshold.data
    #if threshold == None:
    threshold = 0.9
    if "dbtn" in request.form and procopt != "download":
      procopt = 'download'
      return redirect(url_for('range',
        chromosome=form.chromosome.data,
        start=form.start.data,
        end=form.end.data,
        procopt=procopt), code=307)
    if procopt != "download":
      procopt = 'display'
      return redirect(url_for('range',
        chromosome=form.chromosome.data,
        start=form.start.data,
        end=form.end.data,
        procopt=procopt))
  variant_data = []
  if chromosome != None:
    if procopt == 'download':
      logging.info("Make zip data (range)")
      download_list = request.form.getlist("include")
      logging.info("Range - download_list - %s", str(download_list))
      (sample_return_data, snp_return_data, msg) = ds.get_range_data(chromosome, start, end, threshold, download_list)
      zipfilename = "%s_%s_%s_%s" % (chromosome, start, end, "range_results.zip")
      zipFilePath = ds.make_zipfile(sample_return_data, snp_return_data, app.config['UPLOAD_DIR'], zipfilename)
      return send_file(zipFilePath, as_attachment=True)

    (variant_data, msg) = ds.get_variant_data_by_range(chromosome, start, end)
    form.chromosome.data = chromosome
    form.start.data = start
    form.end.data = end
  return render_template('range.html',
    title = 'Range',
    form = form,
    variant_data = variant_data,
    msg = msg,
    db_name = ds.get_db_name())

@app.route('/upload', methods = ['GET', 'POST'])
def upload():
  """
  Display upload form
  """
  form=UploadForm()
  logging.info("Upload [GET]")
  return render_template('upload_file.html', form=form, db_name = ds.get_db_name())

@app.route('/uploaded', methods = ['GET', 'POST'])
@app.route('/uploaded/<filename>/<threshold>', methods = ['GET', 'POST'])
def uploaded_file(filename=None, threshold=0.9):
  """
  Display upload form
  """
  logging.info("uploaded %s", filename)
  form = UploadForm()
  filename = ""
  msg = ""
  if form.validate_on_submit():
    threshold = form.threshold.data
    if threshold == None:
      threshold = 0.9
    filen = form.filename.data
    if filen and allowed_file(filen.filename):
      filename = secure_filename(filen.filename)
      logging.info("uploaded filename = %s (%f)", filename, threshold)
      filen.save(os.path.join(app.config['UPLOAD_DIR'], filename))
      (variant_data, msg) = ds.get_variant_data_for_file(app.config['UPLOAD_DIR'] + "/" + filename, threshold)
      logging.info("variant_data size = %d, msg = %s", len(variant_data), msg)
    dnform=DownloadForm()
    dnform.filename.data=filename
    dnform.filename.readonly=True
    dnform.threshold.data=threshold
    dnform.threshold.readonly=True
    return render_template('uploaded.html',
      title = 'Upload',
      form = dnform,
      variant_data = variant_data,
      msg = msg,
      db_name = ds.get_db_name())
  else:
    logging.info("form validation failed / not run %s, %s", filename, threshold)
  logging.info("uploaded_file - drop through")
  return render_template('upload_file.html', form=form, filename=filename, threshold=threshold)

@app.route('/download', methods = ['POST'])
def download(filename=None, threshold=0.9):
  form=DownloadForm()
  logging.info("download filename(2)= %s", form.filename.data)
  logging.info("download threshold(2) = %s", form.threshold.data)
  threshold = float(form.threshold.data)

  filename = form.filename.data
  download_list = request.form.getlist("include")

  (sample_return_data, snp_return_data, msg) = ds.get_rslist_file_data(app.config['UPLOAD_DIR'] + "/" + filename, threshold, download_list)
  zipfilename = filename + "_file_results.zip"
  zipFilePath = ds.make_zipfile(sample_return_data, snp_return_data, app.config['UPLOAD_DIR'], zipfilename)
  return send_file(zipFilePath, as_attachment=True)
