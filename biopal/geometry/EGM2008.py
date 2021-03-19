# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT


class EGM2008(object):
    def __init__(self, db_path, verbose=False):
        raise NotImplementedError

    def get(self, lat, lon):
        raise NotImplementedError
