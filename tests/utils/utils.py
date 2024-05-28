from io import StringIO
import logging


class CapturingStreamHandler(logging.StreamHandler):
  """Custom StreamHandler to capture logs in memory."""
  def __init__(self):
    super().__init__()
    self.buffer = StringIO()
    self.stream = self.buffer


  def get_logger(self, handler):
    # Get the logger and set the level to capture warnings (adjust if needed)
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        logger.addHandler(handler)
        return logger