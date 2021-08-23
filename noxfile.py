import nox


@nox.session()
def build_sdist(session: nox.Session):
    session.install("build")
    session.run("python", "-m", "build", "--sdist", silent=True)


@nox.session()
def build_wheel(session: nox.Session):
    session.install("build")
    session.run("python", "-m", "build", "--wheel", silent=True)
