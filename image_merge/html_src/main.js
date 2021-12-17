$(document).ready(function () {
    $(window).scroll(function () {
        var scroll = $(window).scrollTop();
        $("#js-hero img").css({
            width: (100 + scroll / 5) + "%"
        })
    })
});