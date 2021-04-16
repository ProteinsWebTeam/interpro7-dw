# -*- coding: utf-8 -*-

from email.message import EmailMessage
from smtplib import SMTP
from typing import Sequence


def notify_curators(host: str, port: int, from_addr: str,
                    to_addr: Sequence[str]):
    msg = EmailMessage()
    msg.set_content("""\
Dear curators,

All steps from the InterPro release procedure requiring data from \
the production database successfully completed.

You may resume integrating signatures, and creating/annotating entries. \
Have fun!

The InterPro Production Team
""")
    msg["Sender"] = from_addr
    msg["To"] = ", ".join(set(to_addr))
    msg["Subject"] = "[curators] Database unfrozen"

    with SMTP(host, port=port) as s:
        s.send_message(msg)
