Integration and Usage
=====================

Integrating into your project
-----------------------------

The fetchscipt can be integrated easily into a project by simply copying the fetchscript directory into the projects root. However, for easy updating it is best to intergrate the fetchscript as a submodule that can be updated from the respository when changes are made.

Usage in a project
------------------

1. Import into your project the classes you need: ``from fetchscript.fetch import FetchDetails``.
2. Initialise an instance on the class: ``f = FetchDetails()``.
3. Call the function to retrieve the information you need: ``entrez_details = f.getDetailsFromEntrez(7157)``.

In case of failure the function will return ``None``.
