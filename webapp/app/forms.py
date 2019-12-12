from flask_wtf import FlaskForm as Form
from wtforms import TextField, FileField, IntegerField, DecimalField
from wtforms.validators import Required, Length, NumberRange, Optional

class SearchForm(Form):
  rs = TextField('rs', validators = [Required()])
  threshold = DecimalField('threshold', validators = [Optional(),NumberRange(min=0.0, max=1.0, message="Must be in range 0.0 to 1.0")])

class SearchRangeForm(Form):
  chr = IntegerField('chr', validators = [Required(),NumberRange(min=1, max=22, message="Must be in range 1 to 22")])
  start = IntegerField('start', validators = [Required(),NumberRange(min=0, message="Must be a +ve integer")])
  end = IntegerField('end', validators = [Required(),NumberRange(min=0, message="Must be a +ve integer")])

class UploadForm(Form):
  filename = FileField('filename', validators = [Required()])
  threshold = DecimalField('threshold', validators = [Optional(),NumberRange(min=0.0, max=1.0, message="Must be in range 0.0 to 1.0")])

class DownloadForm(Form):
  filename = TextField('filename')
  threshold = DecimalField('threshold', validators = [Optional(),NumberRange(min=0.0, max=1.0, message="Must be in range 0.0 to 1.0")])
