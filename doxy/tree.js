var openImg = new Image();
openImg.src = "down.png";
var closedImg = new Image();
closedImg.src = "right.png";

function showBranch(branch){
	var objBranch = document.getElementById(branch).style;
	if(objBranch.display=="block")
		objBranch.display="none";
	else
		objBranch.display="block";
}

function swapFolder(img){
	objImg = document.getElementById(img);
	if(objImg.src.indexOf('right.png')>-1)
		objImg.src = openImg.src;
	else
		objImg.src = closedImg.src;
}

function changeSheets(whichSheet){
    var c = document.Elements.length;
    for(var i=0;i<c;i++)
        var objBranch = document.getElementById(branch).style.display="block";
}

function showAll(){
var arrElements = document.getElementsByClassName("branch");
        for(var i=0; i<arrElements.length; i++){
            arrElements[i].style.display="block";
        }
    }

function closeAll(){
var arrElements = document.getElementsByClassName("branch");
        for(var i=0; i<arrElements.length; i++){
            arrElements[i].style.display="none";
        }
    }

