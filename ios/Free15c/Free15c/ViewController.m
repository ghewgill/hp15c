//
//  ViewController.m
//  Free15c
//
//  Created by Greg Hewgill on 24/11/12.
//  Copyright (c) 2012 Greg Hewgill. All rights reserved.
//

#import "ViewController.h"

@interface ViewController ()

@end

@implementation ViewController {
    UIWebView *core;
    NSMutableDictionary *images;
    UIImageView *calc;
    UIImageView *digit[10];
    UIImageView *decimal[10];
    UIImageView *neg;
    UILabel *user;
    UILabel *f;
    UILabel *g;
    UILabel *trigmode;
    UILabel *complex;
    UILabel *prgm;
    UITapGestureRecognizer *tapper;
}

- (void)viewDidLoad
{
    [super viewDidLoad];
	// Do any additional setup after loading the view, typically from a nib.
    
    NSString *chars = @"0123456789-ABCDEoru";
    images = [[NSMutableDictionary alloc] initWithCapacity:chars.length];
    for (int i = 0; i < chars.length; i++) {
        NSString *c = [chars substringWithRange:NSMakeRange(i, 1)];
        images[c] = [UIImage imageNamed:[NSString stringWithFormat:@"%@.png", c]];
    }
    images[@"."] = [UIImage imageNamed:@"decimal.png"];
    images[@","] = [UIImage imageNamed:@"comma.png"];
    
    calc = [[UIImageView alloc] initWithFrame:self.view.bounds];
    calc.autoresizingMask = UIViewAutoresizingFlexibleWidth | UIViewAutoresizingFlexibleHeight;
    [self.view addSubview:calc];
    //calc.contentMode = UIViewContentModeScaleAspectFit;
    calc.image = [UIImage imageNamed:@"calc.jpg"];
    calc.userInteractionEnabled = YES;
    
    for (int i = 0; i < 10; i++) {
        digit[i] = [[UIImageView alloc] initWithFrame:CGRectMake(65+22 + i * 23, 11+9, 18, 24)];
        digit[i].contentMode = UIViewContentModeTopLeft;
        [self.view addSubview:digit[i]];
        decimal[i] = [[UIImageView alloc] initWithFrame:CGRectMake(65+39 + i * 23, 11+29, 5, 8)];
        decimal[i].contentMode = UIViewContentModeTopLeft;
        [self.view addSubview:decimal[i]];
    }
    
    neg = [[UIImageView alloc] initWithFrame:CGRectMake(65+6, 11+19, 11, 3)];
    neg.hidden = YES;
    neg.contentMode = UIViewContentModeTopLeft;
    neg.image = [UIImage imageNamed:@"neg.png"];
    [self.view addSubview:neg];
    
    UIFont *font = [UIFont systemFontOfSize:9];
    
    user = [[UILabel alloc] initWithFrame:CGRectMake(65+33, 11+37, 26, 8)];
    [self.view addSubview:user];
    user.hidden = YES;
    user.backgroundColor = nil;
    user.opaque = NO;
    user.font = font;
    user.text = @"USER";
    
    f = [[UILabel alloc] initWithFrame:CGRectMake(65+70, 11+37, 4, 8)];
    [self.view addSubview:f];
    f.hidden = YES;
    f.backgroundColor = nil;
    f.opaque = NO;
    f.font = font;
    f.text = @"f";
    
    g = [[UILabel alloc] initWithFrame:CGRectMake(65+84, 11+37, 6, 10)];
    [self.view addSubview:g];
    g.hidden = YES;
    g.backgroundColor = nil;
    g.opaque = NO;
    g.font = font;
    g.text = @"g";
    
    trigmode = [[UILabel alloc] initWithFrame:CGRectMake(65+142, 11+37, 27, 8)];
    [self.view addSubview:trigmode];
    trigmode.hidden = YES;
    trigmode.backgroundColor = nil;
    trigmode.opaque = NO;
    trigmode.font = font;
    trigmode.text = @"";
    
    complex = [[UILabel alloc] initWithFrame:CGRectMake(65+205, 11+37, 7, 8)];
    [self.view addSubview:complex];
    complex.hidden = YES;
    complex.backgroundColor = nil;
    complex.opaque = NO;
    complex.font = font;
    complex.text = @"C";
    
    prgm = [[UILabel alloc] initWithFrame:CGRectMake(65+222, 11+37, 27, 8)];
    [self.view addSubview:prgm];
    prgm.hidden = YES;
    prgm.backgroundColor = nil;
    prgm.opaque = NO;
    prgm.font = font;
    prgm.text = @"PRGM";
    
    tapper = [[UITapGestureRecognizer alloc] initWithTarget:self action:@selector(handleTap:)];
    [calc addGestureRecognizer:tapper];
    
    core = [[UIWebView alloc] init];
    core.delegate = self;
    [core loadRequest:[NSURLRequest requestWithURL:[NSURL fileURLWithPath:[[NSBundle mainBundle] pathForResource:@"core.html" ofType:nil]]]];
}

- (void)didReceiveMemoryWarning
{
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

- (BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)interfaceOrientation
{
    return interfaceOrientation == UIInterfaceOrientationLandscapeLeft
        || interfaceOrientation == UIInterfaceOrientationLandscapeRight;
}

- (BOOL)webView:(UIWebView *)webView shouldStartLoadWithRequest:(NSURLRequest *)request navigationType:(UIWebViewNavigationType)navigationType
{
    //NSLog(@"shouldStartLoadWithRequest %@", request.URL.absoluteString);
    if ([request.URL.scheme isEqualToString:@"callback"]) {
        NSString *decoded = [request.URL.query stringByReplacingPercentEscapesUsingEncoding:NSUTF8StringEncoding];
        NSArray *args = [NSJSONSerialization JSONObjectWithData:[decoded dataUsingEncoding:NSUTF8StringEncoding] options:0 error:nil];
        NSString *method = args.count ? [request.URL.host stringByAppendingString:@":"] : request.URL.host;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Warc-performSelector-leaks"
        [self performSelector:NSSelectorFromString(method) withObject:args];
#pragma clang diagnostic pop
        return NO;
    }
    if ([request.URL.scheme isEqualToString:@"file"]) {
        return YES;
    }
    return NO;
}

- (void)webViewDidFinishLoad:(UIWebView *)webView
{
    [core stringByEvaluatingJavaScriptFromString:@"init()"];
}

- (void)webView:(UIWebView *)webView didFailLoadWithError:(NSError *)error
{
    NSLog(@"didFailLoadWithError %@", error);
}

- (void)clear_digit:(NSArray *)args
{
    int i = [args[0] intValue];
    digit[i].hidden = YES;
}

- (void)clear_digits
{
    for (int i = 0; i < 10; i++) {
        digit[i].hidden = YES;
        decimal[i].hidden = YES;
    }
    neg.hidden = YES;
}

- (void)clear_shift
{
    f.hidden = YES;
    g.hidden = YES;
}

- (void)set_comma:(NSArray *)args
{
    int i = [args[0] intValue];
    decimal[i].image = images[@","];
    decimal[i].hidden = NO;
}

- (void)set_complex:(NSArray *)args
{
    bool on = [args[0] boolValue];
    complex.hidden = !on;
}

- (void)set_decimal:(NSArray *)args
{
    int i = [args[0] intValue];
    decimal[i].image = images[@"."];
    decimal[i].hidden = NO;
}

- (void)set_digit:(NSArray *)args
{
    int i = [args[0] intValue];
    NSString *d = args[1];
    digit[i].image = images[d];
    digit[i].hidden = NO;
}

- (void)set_neg
{
    neg.hidden = NO;
}

- (void)set_prgm:(NSArray *)args
{
    bool on = [args[0] boolValue];
    prgm.hidden = !on;
}

- (void)set_shift:(NSArray *)args
{
    NSString *mode = args[0];
    if ([mode isEqualToString:@"f"]) {
        f.hidden = NO;
    } else if ([mode isEqualToString:@"g"]) {
        g.hidden = NO;
    }
}

- (void)set_trigmode:(NSArray *)args
{
    NSString *mode = args[0];
    if ([mode isEqual:[NSNull null]]) {
        trigmode.hidden = YES;
    } else {
        trigmode.text = mode;
        trigmode.hidden = NO;
    }
}

- (void)set_user:(NSArray *)args
{
    bool on = [args[0] boolValue];
    user.hidden = !on;
}

- (void)handleTap:(UIGestureRecognizer *)sender
{
    CGPoint p = [sender locationInView:calc];
    //NSLog(@"tap %g,%g", p.x, p.y);
    int r = (p.y - 88) * 4 / (300-88);
    int c = (p.x - 0) * 10 / (480-0);
    //NSLog(@"rc %d,%d", r, c);
    if (r < 0 || r > 4 || c < 0 || c > 9) {
        return;
    }
    [core stringByEvaluatingJavaScriptFromString:[NSString stringWithFormat:@"key(KeyTable[%d][%d])", r, c]];
}

@end
