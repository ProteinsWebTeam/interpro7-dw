# -*- coding: utf-8 -*-

from . import ida, relationship


def publish(hosts):
    ida.publish(hosts)
    relationship.publish(hosts)
