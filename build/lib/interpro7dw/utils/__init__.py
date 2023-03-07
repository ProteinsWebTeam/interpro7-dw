import logging
import math


logger = logging.getLogger(__name__)
if not logger.hasHandlers():
    logger.setLevel(logging.INFO)
    _ch = logging.StreamHandler()
    _ch.setFormatter(
        logging.Formatter(
            fmt='%(asctime)s: %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
    )
    logger.addHandler(_ch)


class ProgressTracker:
    def __init__(self, total: int, step: float):
        self.i = 0
        self.total = total
        self.step = self.milestone = math.ceil(total * step)
        self.pc_step = self.pc_milestone = step * 100

    def incr(self):
        self.i += 1

        if self.i == self.milestone:
            logger.info(f"{self.pc_milestone:>10.0f}%")
            self.milestone += self.step

            if self.milestone <= self.total:
                self.pc_milestone += self.pc_step
            else:
                self.milestone = self.total
                self.pc_milestone = 100
