class Logger:
    def __init__(self, quiet = False) -> None:
        self.quiet = quiet
        
    def log(self, message: str) -> None:
        if not self.quiet:
            print(message)