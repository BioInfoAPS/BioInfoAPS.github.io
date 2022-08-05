// Make all tables striped by default.
$("table").addClass("table table-striped");


// Handle foldable challenges and solutions (on click and at start).
$(".solution").click(function(event) {
    var trigger = $(event.target).has(".fold-unfold").length > 0
               || $(event.target).filter(".fold-unfold").length > 0;
    if (trigger) {
        $(">*:not(h2)", this).toggle(400);
        $(">h2>span.fold-unfold", this).toggleClass("glyphicon-collapse-down glyphicon-collapse-up");
        event.stopPropagation();
    }
});
$(".solution").each(function() {
    $(">*:not(h2)", this).toggle();
    var h2 = $("h2:first", this);
    h2.append("<span class='fold-unfold glyphicon glyphicon-collapse-down'></span>");
});


// Handle searches.
// Relies on document having 'meta' element with name 'search-domain'.
function google_search() {
  var query = document.getElementById("google-search").value;
  var domain = $("meta[name=search-domain]").attr("value");
  window.open("https://www.google.com/search?q=" + query + "+site:" + domain);
}

// function to shrink the life cycle bar when scrolling
$(function(){
    $('#life-cycle').data('size','big');
});

$(window).scroll(function(){
    if($(document).scrollTop() > 0)
    {
        if($('#life-cycle').data('size') == 'big')
        {
            $('#life-cycle').data('size','small');
            $('#life-cycle').stop().animate({
                padding: '5px'
            },100);
        }
    }
    else
    {
        if($('#life-cycle').data('size') == 'small')
        {
            $('#life-cycle').data('size','big');
            $('#life-cycle').stop().animate({
                padding: '15px'
            },100);
        }
    }
});

// add copy button to codeblocks

// Copy text as text
function executeCopy(text) {
    var input = document.createElement('textarea');
    document.body.appendChild(input);
    input.value = text;
    input.focus();
    input.select();
    navigator.clipboard.writeText("Hello World!")
    input.remove();
}

// Copy HTML as text (without HTML tags)

document.querySelectorAll(".language-bash").forEach((block) => {
    if (navigator.clipboard) {
        child_code = block.querySelector('pre').innerHTML;
        var doc = new DOMParser().parseFromString(child_code, 'text/html');
        doc.querySelectorAll('.gp').forEach((gp) => {gp.remove()});
        text = doc.body.textContent.replace(/^\s*[\r\n]/gm,"").replace(/^[$>]\s/gm,"");
    
        let button = document.createElement("button");
        button.classList.add('button-cp');
        button.setAttribute('data-cp', text);
        block.appendChild(button);
    }
});
  
document.querySelectorAll('.button-cp').forEach((button) => {
    button.addEventListener('click', () => {
        navigator.clipboard.writeText(button.getAttribute('data-cp'));
        button.classList.add('clicked');
        setTimeout(() => button.classList.remove('clicked'), 1000);
    })
});
