.. highlight:: shell

============
Installation
============

Stable release
--------------

To install teemi, run this command in your terminal:

.. code-block:: console

    $ pip install teemi

This is the preferred method to install teemi, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

From sources
------------

The sources for teemi can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/hiyama341/teemi

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/hiyama341/teemi/tarball/master

Once you have a copy of the source, install `uv <https://docs.astral.sh/uv/getting-started/installation/>`__
and run these commands from the repository root:

.. code-block:: console

    $ uv sync --locked
    $ uv run python

For development, include the development tools and run the tests:

.. code-block:: console

    $ uv sync --locked --extra dev
    $ uv run --locked --extra dev pytest tests

uv manages ``.venv`` and installs Teemi in editable mode. The committed
``uv.lock`` selects dependency versions. Add ``--extra docs`` for documentation
tools or ``--extra learn`` for the H2O notebooks. H2O also requires Java.

Package dependencies are declared in ``pyproject.toml``. Use ``uv add <package>``
to add a runtime dependency and ``uv lock`` after manually editing constraints.
Commit changes to both ``pyproject.toml`` and ``uv.lock``. Use ``uv build`` to
build the source distribution and wheel.

An editable installation with pip is also supported:

.. code-block:: console

    $ pip install -e '.[dev]'

Or directly from PyPI:

.. code-block:: console

    $ pip install 'teemi[dev]'

.. _Github repo: https://github.com/hiyama341/teemi
.. _tarball: https://github.com/hiyama341/teemi/tarball/master
