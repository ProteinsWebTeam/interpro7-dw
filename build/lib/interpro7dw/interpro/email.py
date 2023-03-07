from email.message import EmailMessage
from smtplib import SMTP


def notify_curators(server: str, from_addr: str, to_addrs: str):
    host, port = server.split(":")
    to_addrs = [addr.strip() for addr in to_addrs.split(",")]

    msg = EmailMessage()
    msg.set_content("""\
Dear curators,

All steps from the InterPro release procedure requiring data from \
the production database successfully completed.

You may resume integrating signatures, and creating/annotating entries. \
Have fun!

The InterPro Production Team
""")
    msg["From"] = from_addr
    msg["To"] = ", ".join(set(to_addrs))
    msg["Subject"] = "[curators] Database unfrozen"

    with SMTP(host=host, port=int(port)) as s:
        s.send_message(msg)
