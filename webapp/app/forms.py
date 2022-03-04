from flask_wtf import FlaskForm as Form
from wtforms import StringField, FileField, IntegerField, DecimalField
from wtforms.validators import DataRequired, Length, NumberRange, Optional

class SearchForm(Form):
  rs = StringField('rs', validators = [DataRequired()])
  threshold = DecimalField('threshold', validators = [Optional(),NumberRange(min=0.0, max=1.0, message="Must be in range 0.0 to 1.0")])

class SearchRangeForm(Form):
  chromosome = IntegerField('chromosome', validators = [DataRequired(),NumberRange(min=1, max=22, message="Must be in range 1 to 22")])
  start = IntegerField('start', validators = [DataRequired(),NumberRange(min=0, message="Must be a +ve integer")])
  end = IntegerField('end', validators = [DataRequired(),NumberRange(min=0, message="Must be a +ve integer")])

class UploadForm(Form):
  filename = FileField('filename', validators = [DataRequired()])
  threshold = DecimalField('threshold', validators = [Optional(),NumberRange(min=0.0, max=1.0, message="Must be in range 0.0 to 1.0")])

class DownloadForm(Form):
  filename = StringField('filename')
  threshold = DecimalField('threshold', validators = [Optional(),NumberRange(min=0.0, max=1.0, message="Must be in range 0.0 to 1.0")])
