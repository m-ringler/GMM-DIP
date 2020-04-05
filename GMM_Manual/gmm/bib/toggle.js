/* $Id: toggle.js 133 2007-03-19 17:38:30Z ringler $ */
function toggleAbstract(citekey){
  abs = eval(document.getElementById("abstract_"+citekey));
  if(abs.style.display=="none"){
      abs.style.display="block";
  } else {
    abs.style.display="none";
  }
}