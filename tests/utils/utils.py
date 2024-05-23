from io import StringIO
import logging


class CapturingStreamHandler(logging.StreamHandler):
  """Custom StreamHandler to capture logs in memory."""
  def __init__(self):
    super().__init__()
    self.buffer = StringIO()
    self.stream = self.buffer