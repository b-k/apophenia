s|Outlineheader \([^ ]*\)\(.*\)</p>|<h2><a class="anchor" name="\1"><div class="trigger" onClick="showBranch('\1d');swapFolder('\1f')"><img src="right.png" border="0" id="\1f" alt="pip">\2</div></a></h2><div class="branch" id="\1d">|
s|endofdiv</p>|</div>|
s|ALLBUTTON|<span class="trigger" onClick="showAll();"<a>Expand all </a></span> \| <span class="trigger" onClick="closeAll();"<a>Collapse all </a></span>|
