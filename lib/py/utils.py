from app import app
#------------------------------------------------------
# Allowed-file extension check
# TODO - put this somewhere else!@
#------------------------------------------------------
def allowed_file(filename):
  return '.' in filename and \
    filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

